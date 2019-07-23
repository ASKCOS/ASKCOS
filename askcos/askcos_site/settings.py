"""
Django settings for askcos_site project.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.6/ref/settings/
"""

import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
PROJECT_PATH = os.path.dirname(__file__)

# Get settings from separate celeryconfig.py
from askcos_site.askcos_celery.celeryconfig import *
# Get settings from makeit
import makeit.global_config as gc

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'notsosecret'
# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

ALLOWED_HOSTS = ['0.0.0.0', 'askcos.mit.edu', 'askcos4.mit.edu']
if os.environ.get('CURRENT_HOST'):
    ALLOWED_HOSTS.append(os.environ.get('CURRENT_HOST'))

TEMPLATE_LOADERS = ['django.template.loaders.filesystem.Loader',
 'django.template.loaders.app_directories.Loader']

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            os.path.join(PROJECT_PATH, 'templates'),
        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                # Insert your TEMPLATE_CONTEXT_PROCESSORS here or use this
                # list if you haven't customized them:
                'askcos_site.processors.customization',
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.i18n',
                'django.template.context_processors.media',
                'django.template.context_processors.static',
                'django.template.context_processors.tz',
                'django.template.context_processors.request',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'askcos_site.main',
    'django_extensions',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'askcos_site.urls'
WSGI_APPLICATION = 'askcos_site.wsgi.application'

# LOGIN
LOGIN_URL = '/registration/login'
LOGIN_REDIRECT_URL = '/'
REGISTRATION_OPEN = True
ACCOUNT_ACTIVATION_DAYS=7
REGISTRATION_SALT='saltystring'

# Registration
REGISTRATION_SUPPLEMENT_CLASS = None
ACCOUNT_ACTIVATION_DAYS = 7
# Uses postfix server running at localhost:25
#EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'
#EMAIL_HOST = 'localhost'
#    EMAIL_PORT = 25
#EMAIL_HOST_USER = ''
#EMAIL_HOST_PASSWORD = ''
#EMAIL_USE_TLS = False
DEFAULT_FROM_EMAIL = 'no-reply@askcos4.mit.edu'

# Where are user settings / banlists / etc. saved?
# NOTE: we recommend relocating the db to an ssd for speed
DATABASES = {'default': {
    'ENGINE': 'django.db.backends.sqlite3',
    'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
}}


LANGUAGE_CODE = 'en-us'
TIME_ZONE = 'UTC'
USE_I18N = True
USE_L10N = True
USE_TZ = True

# Static files (CSS, JavaScript, Images)
STATIC_ROOT = os.path.join(PROJECT_PATH, 'static/')
STATIC_URL = '/static/'
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',    
) 
STATICFILES_DIRS = (
    os.path.join(STATIC_ROOT, 'css'),
    os.path.join(STATIC_ROOT, 'js')
)

# Media files
MEDIA_ROOT = os.path.join(PROJECT_PATH, 'media/')
MEDIA_URL = '/media/'

################################################################################
# Define databases to replicate Make-It settings
################################################################################

RETRO_TRANSFORMS = gc.RETRO_TRANSFORMS
RETRO_TRANSFORMS_CHIRAL = gc.RETRO_TRANSFORMS_CHIRAL

# TODO: deprecate in favor of template-free forward prediction only
SYNTH_TRANSFORMS = gc.SYNTH_TRANSFORMS

INSTANCES = gc.INSTANCES
REACTIONS = gc.REACTIONS
CHEMICALS = gc.CHEMICALS

# Back-up databases
INSTANCES_OLD = {'database': 'reaxys', 'collection': 'instances'}
REACTIONS_OLD = {'database': 'reaxys', 'collection': 'reactions'}
CHEMICALS_OLD = {'database': 'reaxys', 'collection': 'chemicals'}

BUYABLES = gc.BUYABLES
SOLVENTS = gc.SOLVENTS

# For searching "old" templates in the MongoDB (not locally saved)
TEMPLATE_BACKUPS = [
    ('reaxys_v2', 'transforms_retro_v8'),
    ('reaxys_v2', 'transforms_retro_v7'),
    ('reaxys', 'transforms_retro_v6'),
    ('reaxys', 'transforms_retro_v5'),
    ('reaxys', 'transforms_retro_v4'),
    ('reaxys', 'transforms_retro_v3'),
]

################################################################################
# Define local file storage locations
################################################################################

### Very important - where to look for local versions of files instead of relying on DB connections
LOCAL_STORAGE = {}
LOCAL_STORAGE['user_saves'] = os.path.join(gc.data_path, 'user_saves')


