def get_retro_transformer():
    '''
    Loads the retro transformer
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
    # v3 uses rdchiral, for chiral transformations
    import makeit.webapp.transformer_v3 as transformer
    RetroTransformer = transformer.Transformer()
    mincount_retro = settings.RETRO_TRANSFORMS_CHIRAL['mincount']
    mincount_retro_chiral = settings.RETRO_TRANSFORMS_CHIRAL['mincount_chiral']
    RetroTransformer.load(RETRO_DB, mincount=mincount_retro, get_retro=True, 
        get_synth=False, mincount_chiral=mincount_retro_chiral)
    print('Chiral retro transformer loaded {} retro templates'.format(RetroTransformer.num_templates))

    # Also - get the Pricer and load it into the RetroTransformer object.
    db = db_client[settings.BUYABLES['database']]
    BUYABLE_DB = db[settings.BUYABLES['collection']]
    db = db_client[settings.CHEMICALS['database']]
    CHEMICAL_DB = db[settings.CHEMICALS['collection']]
    print('Loading prices...')
    import makeit.retro.pricer as pricer
    RetroTransformer.Pricer = pricer.Pricer()
    RetroTransformer.Pricer.load(CHEMICAL_DB, BUYABLE_DB)
    print('Loaded known prices')

    return RetroTransformer