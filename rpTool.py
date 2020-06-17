import copy
import logging


## Class to add the cofactors to a monocomponent reaction to construct the full reaction
#
#
class rpCofactors:
    ## Init method
    # Here we want to seperate what is the use input and what is parsed by the cache to make sure that
    # everything is not hadled by a single
    #
    # @param rpReader input reader object with the parsed user input and cache files required
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpCofactors')
        ##### stuff to load from cache #####
        self.deprecatedCID_cid = {}
        self.deprecatedRID_rid = {}
        self.cid_strc = None
        self.rr_full_reactions = None
        self.cid_xref = None
        self.rr_reactions = None


    ################################################################
    ######################## PRIVATE FUNCTIONS #####################
    ################################################################


    ## Function to create a dictionnary of old to new reaction id's
    #
    def _checkRIDdeprecated(self, rid):
        try:
            return self.deprecatedRID_rid[rid]
        except KeyError:
            return rid

    ## Function to create a dictionnary of old to new chemical id's
    #
    #  Generate a one-to-one dictionnary of old id's to new ones. Private function
    #
    def _checkCIDdeprecated(self, cid):
        try:
            return self.deprecatedCID_cid[cid]
        except KeyError:
            return cid


    ################################################################
    ######################### PUBLIC FUNCTIONS #####################
    ################################################################


    ## Given a dictionnary describing a monocomponent reaction, add the cofactors by comparing it with the original reaction
    #
    # @param step Dictionnary describing the reaction
    # @param reac_side String 'right' or 'left' describing the direction of the monocomponent reaction compared with the original reaction
    # @param rr_reac Dictionnary describing the monocomponent reaction from RetroRules
    # @param f_reac Dictionnary describing the full original reaction
    # @param pathway_cmp Dictionnary used to retreive the public ID of the intermediate compounds. Resets for each individual pathway
    #
    def completeReac(self, step, rr_reac, full_reac, mono_side, rr_string, pathway_cmp):
        if mono_side:
            ## add the unknown species to pathway_cmp for the next steps
            rr_mono_cmp = list(rr_reac.keys())
            step_mono_cmp = list(step.keys())
            if (len(rr_mono_cmp)==1 and len(step_mono_cmp)==1):
                #this is purposely overwitten since the main cmp between reactions can change
                pathway_cmp[step_mono_cmp[0]] = rr_mono_cmp[0]
            else:
                self.logger.warning('There should be only one compound on the left for monocomponent reaction: rr_mono_cmp: '+str(rr_mono_cmp)+' step_mono_cmp: '+str(step_mono_cmp))
                return False
        ## add the side species
        for toAdd in full_reac.keys()-rr_reac.keys():
            step.update({toAdd: full_reac[toAdd]})
            ### update the reaction rule string
            try:
                smi = self.cid_strc[toAdd]['smiles']
                if not smi==None:
                    for sto_add in range(int(full_reac[toAdd])):
                        rr_string += '.'+str(smi)
            except KeyError:
                self.logger.warning('Cannot find smiles structure for '+str(toAdd))
        ## Update the the stochio
        for step_spe in step:
            if step_spe in full_reac:
                if not step[step_spe]==full_reac[step_spe]:
                    stochio_diff = full_reac[step_spe]-step[step_spe]
                    step[step_spe] = full_reac[step_spe]
                    if stochio_diff<0:
                        self.logger.warning('full_reac stochio should never be smaller than step')
                        continue
                    for i in range(stochio_diff):
                        ### update the reaction rule string
                        try:
                            smi = self.cid_strc[step_spe]['smiles']
                            if not smi==None:
                                rr_string += '.'+str(smi)
                        except KeyError:
                            self.logger.warning('Cannot find smiles structure for '+str(toAdd))
            elif step_spe in pathway_cmp:
                if pathway_cmp[step_spe] in full_reac:
                    if not step[step_spe]==full_reac[pathway_cmp[step_spe]]:
                        step[step_spe] = full_reac[pathway_cmp[step_spe]]
            #Its fine if the stochio is not updated, better than ignoring a whole pathway
                #else:
                #    self.logger.warning('Cannot find '+str(step_spe)+' in full reaction')
                #    return False
            #else:
            #    self.logger.warning('Cannot find '+str(step_spe)+' in pathway_cmp')
            #    return False
        return True, rr_string


    ## Add the cofactors to monocomponent reactions
    #
    # @param step Step in a pathway
    # @param pathway_cmp Dictionnary of intermediate compounds with their public ID's
    # @return Boolean determine if the step is to be added
    def addCofactors_step(self, step, pathway_cmp):
        reac_smiles_left = step['reaction_rule'].split('>>')[0]
        reac_smiles_right = step['reaction_rule'].split('>>')[1]
        if self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']==-1:
            isSuccess, reac_smiles_left = self.completeReac(step['right'],
                    self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['left'],
                    self.rr_full_reactions[step['rule_ori_reac']]['right'],
                    True,
                    reac_smiles_left,
                    pathway_cmp)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
            isSuccess, reac_smiles_right = self.completeReac(step['left'],
                    self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['right'],
                    self.rr_full_reactions[step['rule_ori_reac']]['left'],
                    False,
                    reac_smiles_right,
                    pathway_cmp)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
        elif self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']==1:
            isSuccess, reac_smiles_left = self.completeReac(step['right'],
                    self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['left'],
                    self.rr_full_reactions[step['rule_ori_reac']]['left'],
                    True,
                    reac_smiles_left,
                    pathway_cmp)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
            isSuccess, reac_smiles_right = self.completeReac(step['left'],
                    self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['right'],
                    self.rr_full_reactions[step['rule_ori_reac']]['right'],
                    False,
                    reac_smiles_right,
                    pathway_cmp)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
        else:
            self.logger.error('Relative direction can only be 1 or -1: '+str(self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']))
            return False
        step['reaction_rule'] = reac_smiles_left+'>>'+reac_smiles_right
        return True


    ## Function to reconstruct the heterologous pathway
    #
    #  Read each pathway information and RetroRules information to construct heterologous pathways and add the cofactors
    #
    #  @param self Object pointer
    #  @param rpsbml rpSBML object with a single model
    #  @return Boolean if True then you keep that model for the next step, if not then ignore it
    def addCofactors(self, rpsbml, compartment_id='MNXC3', pathway_id='rp_pathway'):
        #This keeps the IDs conversions to the pathway
        pathway_cmp = {}
        rp_path = rpsbml.outPathsDict(pathway_id)
        ori_rp_path = copy.deepcopy(rp_path)
        #We reverse the loop to ID the intermediate CMP to their original ones
        for stepNum in sorted(list(rp_path), reverse=True):
        #for stepNum in sorted(list(rp_path)):
            if self.addCofactors_step(rp_path[stepNum], pathway_cmp):
                ###add the new cofactors to the SBML
                #remove the original species from the monocomponent reaction
                reactants = set(set(rp_path[stepNum]['left'].keys())-set(ori_rp_path[stepNum]['left'].keys()))
                products = set(set(rp_path[stepNum]['right'].keys())-set(ori_rp_path[stepNum]['right'].keys()))
                for species in reactants|products:
                    #check to make sure that they do not yet exist and if not create a new one
                    #TODO, replace the species with an existing one if it is contained in the MIRIAM annotations
                    tmp_species = self._checkCIDdeprecated(species)
                    #neeed to test all the MIRIAM species comparison
                    if not rpsbml.speciesExists(tmp_species, compartment_id):
                        xref = {}
                        inchi = None
                        inchikey = None
                        smiles = None
                        chemName = None
                        try:
                            xref = self.cid_xref[tmp_species]
                        except KeyError:
                            try:
                                xref = self.cid_xref[tmp_species]
                            except KeyError:
                                #TODO: although there should not be any
                                #intermediate species here consider
                                #removing this warning
                                self.logger.warning('Cannot find the xref for this species: '+str(tmp_species))
                                pass
                        try:
                            inchi = self.cid_strc[tmp_species]['inchi']
                        except KeyError:
                            try:
                                inchi = self.cid_strc[tmp_species]['inchi']
                            except KeyError:
                                self.logger.warning('Cannot find the inchi for this species: '+str(tmp_species))
                                pass
                        try:
                            inchikey = self.cid_strc[tmp_species]['inchikey']
                        except KeyError:
                            self.logger.warning('Cannot find the inchikey for this species: '+str(tmp_species))
                            pass
                        try:
                            smiles = self.cid_strc[tmp_species]['smiles']
                        except KeyError:
                            self.logger.warning('Cannot find the smiles for this species: '+str(tmp_species))
                            pass
                        #add the new species to rpsbml
                        try:
                            chemName = self.cid_strc[tmp_species]['name']
                        except KeyError:
                            self.logger.warning('Cannot find the name for this species: '+str(tmp_species))
                            pass
                        rpsbml.createSpecies(tmp_species,
                                             compartment_id,
                                             chemName,
                                             xref,
                                             inchi,
                                             inchikey,
                                             smiles)
                #add the new species to the RP reactions
                reac = rpsbml.model.getReaction(rp_path[stepNum]['reaction_id'])
                for pro in products:
                    prod = reac.createProduct()
                    prod.setSpecies(str(self._checkCIDdeprecated(pro))+'__64__'+str(compartment_id))
                    prod.setConstant(True)
                    prod.setStoichiometry(rp_path[stepNum]['right'][pro])
                for sub in reactants:
                    subs = reac.createReactant()
                    subs.setSpecies(str(self._checkCIDdeprecated(sub))+'__64__'+str(compartment_id))
                    subs.setConstant(True)
                    subs.setStoichiometry(rp_path[stepNum]['left'][sub])
                #replace the reaction rule with new one
                rpsbml.addUpdateBRSynth(reac, 'smiles', rp_path[stepNum]['reaction_rule'], None, True)
            else:
                #if the cofactors cannot be found delete it from the list
                self.logger.warning('Cannot find cofactors... skipping')
                return False
        return True
