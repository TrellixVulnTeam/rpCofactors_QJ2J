import copy
import logging
import time
import json


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
        self.cid_name = None
        self.inchikey_cid = None
        self.rr_reactions = None
        ##### pubchem search ###
        self.pubchem_inchi = {}
        self.pubchem_inchikey = {}
        self.pubchem_smiles = {} 
        self.pubchem_min_count = 0
        self.pubchem_min_start = 0.0



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

    ## Function to create return the uniform compound ID
    #
    # @param cid String The identifier for a given compounf
    #
    def _checkCIDdeprecated(self, cid):
        try:
            return self.deprecatedCID_cid[cid]
        except KeyError:
            return cid


    ## Function that waits for the time limit of the pubchem web service if exceeded
    #
    #
    def _pubChemLimit(self):
        if self.pubchem_min_start==0.0:
            self.pubchem_min_start = time.time()
        #self.pubchem_sec_count += 1
        self.pubchem_min_count += 1
        #### requests per minute ####
        if self.pubchem_min_count>=500 and time.time()-self.pubchem_min_start<=60.0:
            self.logger.warning('Reached 500 requests per minute for pubchem... waiting a minute')
            time.sleep(60.0)
            self.pubchem_min_start = time.time()
            self.pubchem_min_count = 0
        elif time.time()-self.pubchem_min_start>60.0:
            self.pubchem_min_start = time.time()
            self.pubchem_min_count = 0

    ## Try to retreive the xref from an inchi structure using pubchem
    #
    #
    '''
    No more than 5 requests per second.
    No more than 400 requests per minute.
    No longer than 300 second running time per minute.
    Requests exceeding limits are rejected (HTTP 503 error)
    '''
    def _pubchemStrctSearch(self, strct, itype='inchi'):
        self._pubChemLimit()
        try:
            r = requests.post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'+str(itype)+'/xrefs/SBURL/JSON', data={itype: strct})
            res_list = r.json()
        except json.decoder.JSONDecodeError:
            self.logger.warning('JSON decode error')
            return {}
        try:
            res_list = res_list['InformationList']['Information']
        except KeyError:
            self.logger.warning('pubchem JSON keyerror: '+str(res_list))
            return {}
        xref = {}
        if len(res_list)==1:
            self._pubChemLimit()
            try:
                prop = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(res_list[0]['CID'])+'/property/IUPACName,InChI,InChIKey,CanonicalSMILES/JSON')
                prop_list = prop.json()
            except json.decoder.JSONDecodeError:
                self.logger.warning('JSON decode error')
                return {}
            try:
                name = prop_list['PropertyTable']['Properties'][0]['IUPACName']
                inchi = prop_list['PropertyTable']['Properties'][0]['InChI']
                inchikey = prop_list['PropertyTable']['Properties'][0]['InChIKey']
                smiles = prop_list['PropertyTable']['Properties'][0]['CanonicalSMILES']
            except KeyError:
                self.logger.warning('pubchem JSON keyerror: '+str(prop_list))
                return {}
            #TODO: need to determine how long cobra cannot handle this
            #TODO: determine if names that are too long is the problem and if not remove this part
            if len(name)>30:
                self._pubChemLimit()
                try:
                    syn = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(res_list[0]['CID'])+'/synonyms/JSON')
                    syn_lst = syn.json()
                except json.decoder.JSONDecodeError:
                    self.logger.warning('pubchem JSON decode error')
                    return {}
                try:
                    syn_lst = syn_lst['InformationList']['Information'][0]['Synonym']
                    syn_lst = [x for x in syn_lst if not 'CHEBI' in x and not x.isupper()]
                    name = syn_lst[0] #need a better way instead of just the firs tone
                except KeyError:
                    self.logger.warning('pubchem JSON keyerror: '+str(syn.json()))
                    return {}
                except IndexError:
                    name = ''
            xref['pubchem'] = [str(res_list[0]['CID'])]
            for url in res_list[0]['SBURL']:
                if 'https://biocyc.org/compound?orgid=META&id=' in url:
                    if 'biocyc' not in xref:
                        xref['biocyc'] = []
                    xref['biocyc'].append(url.replace('https://biocyc.org/compound?orgid=META&id=', ''))
                if 'http://www.hmdb.ca/metabolites/' in url:
                    if 'hmdb' not in xref:
                        xref['hmdb'] = []
                    xref['hmdb'].append(url.replace('http://www.hmdb.ca/metabolites/', ''))
                if 'http://www.genome.jp/dbget-bin/www_bget?cpd:' in url:
                    if 'kegg_c' not in xref:
                        xref['kegg_c'] = []
                    xref['kegg_c'].append(url.replace('http://www.genome.jp/dbget-bin/www_bget?cpd:', ''))
                if 'http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' in url:
                    if 'chebi' not in xref:
                        xref['chebi'] = []
                    xref['chebi'].append(url.replace('http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:', ''))
        elif len(res_list)==0:
            self.logger.warning('Could not find results for: '+str(strct))
            return {}
        else:
            self.logger.warning('There are more than one result for '+str(strct)+'... Ignoring')
            return {}
        return {'name': name, 'inchi': inchi, 'inchikey': inchikey, 'smiles': smiles, 'xref': xref}


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
    def addCofactors(self, rpsbml, compartment_id='MNXC3', pathway_id='rp_pathway', pubchem_search=False):
        #This keeps the IDs conversions to the pathway
        pathway_cmp = {}
        rpsbml_json = rpsbml.genJSON(pathway_id)
        rp_path = rpsbml.outPathsDict(pathway_id)
        ori_rp_path = copy.deepcopy(rp_path)
        self.logger.debug(rpsbml_json)
        spe_conv = {}
        #We reverse the loop to ID the intermediate CMP to their original ones
        for stepNum in sorted(list(rp_path), reverse=True):
        #for stepNum in sorted(list(rp_path)):
            if self.addCofactors_step(rp_path[stepNum], pathway_cmp):
                ###add the new cofactors to the SBML
                #remove the original species from the monocomponent reaction
                reactants = set(set(rp_path[stepNum]['left'].keys())-set(ori_rp_path[stepNum]['left'].keys()))
                products = set(set(rp_path[stepNum]['right'].keys())-set(ori_rp_path[stepNum]['right'].keys()))
                for species in reactants|products:
                    tmp_species = self._checkCIDdeprecated(species)
                    self.logger.debug('----------- '+str(tmp_species)+' -------------')
                    self.logger.debug('spe_conv: '+str(spe_conv))
                    if not tmp_species in spe_conv and not rpsbml.speciesExists(tmp_species, compartment_id):
                        #check to make sure that they do not yet exist and if not create a new one
                        #TODO, replace the species with an existing one if it is contained in the MIRIAM annotations
                        ##### 1) gather the species information
                        xref = {}
                        inchi = None
                        inchikey = None
                        smiles = None
                        chem_name = None
                        try:
                            inchi = self.cid_strc[tmp_species]['inchi']
                        except KeyError:
                            self.logger.warning('Cannot find the inchi for this species: '+str(tmp_species))
                        try:
                            inchikey = self.cid_strc[tmp_species]['inchikey']
                            self.logger.debug('Found the inchikey: '+str(inchikey))
                            #### TODO: find a better way to check if two species are the same ####
                            isfound = False
                            for rpsbml_species in rpsbml_json['species']:
                                #TODO add a comparison by xref as well
                                #self.logger.debug(str(rpsbml_json['species'][rpsbml_species]['brsynth']['inchikey'])+' <--> '+str(inchikey))
                                if str(rpsbml_json['species'][rpsbml_species]['brsynth']['inchikey'])==str(inchikey):
                                    spe_conv[tmp_species] = rpsbml_species
                                    self.logger.debug('The species '+str(tmp_species)+' is the same as '+str(rpsbml_species))
                                    isfound = True
                                    break
                            if isfound:
                                continue
                        except KeyError:
                            self.logger.warning('Cannot find the inchikey for this species: '+str(tmp_species))
                        try:
                            smiles = self.cid_strc[tmp_species]['smiles']
                        except KeyError:
                            self.logger.warning('Cannot find the smiles for this species: '+str(tmp_species))
                        try:
                            xref = self.cid_xref[tmp_species]
                        except KeyError:
                            try:
                                xref = self.cid_xref[tmp_species]
                            except KeyError:
                                #if you cannot find using cid, try to retreive it using its inchikey
                                try:
                                    if inchikey:
                                        # WARNING here we use MNX since as of now, there are only MNX data that is parsed correctly
                                        tmp_cids = [i for i in self.inchikey_cid[inchikey] if i[:3]=='MNX']
                                        #TODO: handle multiple matches. For now we assume that multiple MNX means that there are deprecated versions of the tool
                                        if tmp_cids:
                                            xref = self.cid_xref[self._checkCIDdeprecated(tmp_cids[0])]
                                except KeyError:
                                    self.logger.warning('Cannot find the xref for this species: '+str(tmp_species))
                                    xref = {}
                        #### Common Name ####
                        try:
                            chem_name = self.cid_name[self._checkCIDdeprecated(tmp_species)]
                        except KeyError:
                            #if you cannot find using cid, try to retreive it using its inchikey
                            try:
                                if inchikey:
                                    tmp_cids = [i for i in self.inchikey_cid[inchikey] if i[:3]=='MNX']
                                    if tmp_cids:
                                        chem_name = self.cid_name[self._checkCIDdeprecated(tmp_cids[0])]
                            except KeyError:
                                self.logger.warning('Cannot find the name for this species: '+str(tmp_species))
                        ###### Try to recover information using the structures and pubchem REST requests ####
                        pubchem_inchi = None
                        pubchem_inchikey = None
                        pubchem_smiles = None
                        pubchem_xref = {}
                        #inchi
                        try:
                            if not xref and pubchem_search:
                                try:
                                    pubchem_inchi = self.pubchem_inchi[inchi]['inchi']
                                    pubchem_inchikey = self.pubchem_inchi[inchi]['inchikey']
                                    pubchem_smiles = self.pubchem_inchi[inchi]['smiles']
                                    pubchem_xref = self.pubchem_inchi[inchi]['xref'] 
                                except KeyError:
                                    #if not self.pubchem_inchi[inchi]=={}:
                                    pubres = self._pubchemStrctSearch(inchi, 'inchi')
                                    if not chem_name:
                                        chem_name = pubres['name']
                                    if 'chebi' in pubres['xref']:
                                        try:
                                            xref = self.cid_xref[self.chebi_cid[pubres['xref']['chebi'][0]]]
                                        except KeyError:
                                            pass
                                    if not pubchem_xref:
                                        pubchem_xref = pubres['xref']
                                    if not pubchem_inchikey:
                                        pubchem_inchikey = pubres['inchikey']
                                    if not pubchem_smiles:
                                        pubchem_smiles = pubres['smiles']
                        except KeyError:
                            self.logger.warning('Bad results from pubchem results')
                            self.pubchem_inchi[inchi] = {}
                            pass
                        #inchikey
                        try:
                            if not xref and pubchem_search:
                                try:
                                    pubchem_inchi = self.pubchem_inchikey[inchikey]['inchi']
                                    pubchem_inchikey = self.pubchem_inchikey[inchikey]['inchikey']
                                    pubchem_smiles = self.pubchem_inchikey[inchikey]['smiles']
                                    pubchem_xref = self.pubchem_inchikey[inchikey]['xref'] 
                                except KeyError:
                                    #if not self.pubchem_inchikey[inchikey]=={}:
                                    pubres = self._pubchemStrctSearch(inchikey, 'inchikey')
                                    if not chem_name:
                                        chem_name = pubres['name']
                                    if 'chebi' in pubres['xref']:
                                        try:
                                            xref = self.cid_xref[self.chebi_cid[pubres['xref']['chebi'][0]]]
                                        except KeyError:
                                            pass
                                    if not pubchem_xref:
                                        pubchem_xref = pubres['xref']
                                    if not pubchem_inchi:
                                        pubchem_inchi = pubres['inchi']
                                    if not pubchem_smiles:
                                        pubchem_smiles = pubres['smiles']
                        except KeyError:
                            self.logger.warning('Bad results from pubchem results')
                            self.pubchem_inchikey[inchikey] = {}
                            pass
                        #smiles
                        try:
                            if not xref and pubchem_search:
                                try:
                                    pubchem_inchi = self.pubchem_smiles[smiles]['inchi']
                                    pubchem_inchikey = self.pubchem_smiles[smiles]['inchikey']
                                    pubchem_smiles = self.pubchem_smiles[smiles]['smiles']
                                    pubchem_xref = self.pubchem_smiles[smiles]['xref'] 
                                except KeyError:
                                    #if not self.pubchem_smiles[smiles]=={}:
                                    pubres = self._pubchemStrctSearch(smiles, 'smiles')
                                    if not chem_name:
                                        chem_name = pubres['name']
                                    if 'chebi' in pubres['xref']:
                                        try:
                                            xref = self.cid_xref[self.chebi_cid[pubres['xref']['chebi'][0]]]
                                        except KeyError:
                                            pass
                                    if not pubchem_xref:
                                        pubchem_xref = pubres['xref']
                                    if not pubchem_inchi:
                                        pubchem_inchi = pubres['inchi']
                                    if not pubchem_inchikey:
                                        pubchem_inchikey = pubres['inchikey']
                        except KeyError:
                            self.pubchem_smiles[smiles] = {}
                            self.logger.warning('Bad results from pubchem results')
                            pass
                        if not inchi:
                            inchi = pubchem_inchi
                        if not inchikey:
                            inchikey = pubchem_inchikey
                        if not smiles:
                            smiles = pubchem_smiles
                        if pubchem_inchi:
                            self.pubchem_inchi[pubchem_inchi] = {'inchi': pubchem_inchi, 'smiles': pubchem_smiles, 'inchikey': pubchem_inchikey, 'xref': pubchem_xref}
                        if pubchem_inchikey:
                            self.pubchem_inchikey[pubchem_inchikey] = {'inchi': pubchem_inchi, 'smiles': pubchem_smiles, 'inchikey': pubchem_inchikey, 'xref': pubchem_xref}
                        if pubchem_smiles:
                            self.pubchem_smiles[pubchem_smiles] = {'inchi': pubchem_inchi, 'smiles': pubchem_smiles, 'inchikey': pubchem_inchikey, 'xref': pubchem_xref}
                        if not xref:
                            xref = pubchem_xref
                        #pass the information to create the species
                        if chem_name:
                            chem_name = chem_name.replace("'", "")
                        #add the new species to rpsbml
                        #### TODO: find a better way to check if two species are the same ####
                        for rpsbml_species in rpsbml_json['species']:
                            #TODO add a comparison by xref as well
                            ##self.logger.debug(str(rpsbml_json['species'][rpsbml_species]['brsynth']['inchikey'])+' <--> '+str(inchikey))
                            isfound = False
                            if str(rpsbml_json['species'][rpsbml_species]['brsynth']['inchikey'])==str(inchikey):
                                spe_conv[tmp_species] = rpsbml_species
                                self.logger.debug('The species '+str(tmp_species)+' is the same as '+str(rpsbml_species))
                                isfound = True
                                break
                            if isfound:
                                continue
                        self.logger.debug('Creating species: '+str(tmp_species))
                        rpsbml.createSpecies(tmp_species,
                                             compartment_id,
                                             chem_name,
                                             xref,
                                             inchi,
                                             inchikey,
                                             smiles) 
                #add the new species to the RP reactions
                reac = rpsbml.model.getReaction(rp_path[stepNum]['reaction_id'])
                for pro in products:
                    prod = reac.createProduct()
                    if self._checkCIDdeprecated(pro) in spe_conv:
                        toadd = spe_conv[self._checkCIDdeprecated(pro)]
                    else:
                        toadd = str(self._checkCIDdeprecated(pro))+'__64__'+str(compartment_id)
                    #prod.setSpecies(str(self._checkCIDdeprecated(pro))+'__64__'+str(compartment_id))
                    prod.setSpecies(toadd)
                    prod.setConstant(True)
                    prod.setStoichiometry(rp_path[stepNum]['right'][pro])
                for sub in reactants:
                    subs = reac.createReactant()
                    if self._checkCIDdeprecated(sub) in spe_conv:
                        toadd = spe_conv[self._checkCIDdeprecated(sub)]
                    else:
                        toadd = str(self._checkCIDdeprecated(pro))+'__64__'+str(compartment_id)
                    #prod.setSpecies(str(self._checkCIDdeprecated(sub))+'__64__'+str(compartment_id))
                    prod.setSpecies(toadd)
                    subs.setConstant(True)
                    subs.setStoichiometry(rp_path[stepNum]['left'][sub])
                #replace the reaction rule with new one
                rpsbml.addUpdateBRSynth(reac, 'smiles', rp_path[stepNum]['reaction_rule'], None, True)
            else:
                #if the cofactors cannot be found delete it from the list
                self.logger.warning('Cannot find cofactors... skipping')
                return False
        return True
