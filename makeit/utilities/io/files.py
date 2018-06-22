import os
import makeit.global_config as gc 

def make_directory(dir_name):
    path = os.path.join(os.getcwd(), dir_name)
    if not os.path.isdir(path):
        os.mkdir(path)
    return path

################################################################################
# Where are local files stored?
################################################################################

def get_retrotransformer_achiral_path(dbname, collname, mincount_retro):
    return os.path.join(settings.LOCAL_STORAGE['dir'], 
        'retrotransformer_achiral_using_%s-%s_mincount%i.pkl' % (dbname, collname, mincount_retro))

def get_retrotransformer_chiral_path(dbname, collname, mincount_retro, mincount_retro_chiral):
    return os.path.join(settings.LOCAL_STORAGE['dir'], 
        'retrotransformer_chiral_using_%s-%s_mincount%i_mincountchiral%i.pkl' % (dbname, collname, mincount_retro, mincount_retro_chiral))

def get_synthtransformer_path(dbname, collname, mincount):
    return os.path.join(settings.LOCAL_STORAGE['dir'], 
        'synthtransformer_using_%s-%s_mincount%i.pkl' % (dbname, collname, mincount))

def get_pricer_path(chem_dbname, chem_collname, buyable_dbname, buyable_collname):
    return os.path.join(settings.LOCAL_STORAGE['dir'], 
        'pricer_using_%s-%s_and_%s-%s.pkl' % (chem_dbname, chem_collname, buyable_dbname, buyable_collname))
