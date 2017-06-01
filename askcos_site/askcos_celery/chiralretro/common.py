def get_retro_transformer():
    '''
    Loads the retro transformer that uses Stereofix
    '''
    print('Loading chiral retro transofmer')

    # Get Django settings
    from django.conf import settings

    # Database
    from database import db_client
    db = db_client[settings.INSTANCES['database']]
    RETRO_DB = db[settings.RETRO_TRANSFORMS_CHIRAL['collection']]

    # Setting logging low
    from rdkit import RDLogger
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    # Import retro transformer class and load
    # v2 uses Stereofix, which is necessary for chiral transformations
    import makeit.webapp.transformer_v2 as transformer
    RetroTransformer = transformer.Transformer()
    mincount_retro = settings.RETRO_TRANSFORMS_CHIRAL['mincount']
    mincount_retro_chiral = settings.RETRO_TRANSFORMS_CHIRAL['mincount_chiral']
    RetroTransformer.load(RETRO_DB, mincount=mincount_retro, get_retro=True, 
        get_synth=False, mincount_chiral=mincount_retro_chiral)
    print('Chiral retro transformer loaded {} retro templates'.format(RetroTransformer.num_templates))
    return RetroTransformer