import rdkit.Chem as Chem 
import rdkit.Chem.AllChem as AllChem 

def load_templates(SYNTH_DB, mincount=0, countsonly=False):
    print('Loading synthetic templates for forward predictor')
    templates = []
    for doc in SYNTH_DB.find({'count': {'$gte': mincount}}, ['_id', 'count', 'reaction_smarts']):
        if 'reaction_smarts' not in doc: continue
        reaction_smarts = str(doc['reaction_smarts'])
        if not reaction_smarts: continue
        template = {
            'reaction_smarts':      reaction_smarts,
            'count':                doc['count'] if 'count' in doc else 0,
            '_id':                  doc['_id'] if '_id' in doc else -1,
        }
        try:
            rxn_f = AllChem.ReactionFromSmarts('(' + reaction_smarts.replace('>>', ')>>(') + ')')
            if rxn_f.Validate()[1] != 0:
                print('Could not validate {}'.format(reaction_smarts))
                continue
            template['rxn_f'] = rxn_f
        except Exception as e:
            print('Couldnt load forward: {}: {}'.format(reaction_smarts, e))
            continue
        templates.append(template)

    num_templates = len(templates)
    templates = sorted(templates, key=lambda z: z['count'], reverse=True)
    print('Loaded {} templates'.format(num_templates))

    if countsonly:
        return [x['count'] for x in templates]
    return templates