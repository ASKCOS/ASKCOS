import os
import sys


os.environ['DJANGO_SETTINGS_MODULE'] = 'askcos_site.settings'


from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()

# Load libcairo
from ctypes import cdll
cdll.LoadLibrary('libcairo.so.2')
import cairo 
cairo_dir = dir(cairo)
if len(cairo_dir) < 20:
    print('######################################')
    print('## CAIRO INSTALLATION LIKELY BROKEN ##')
    print('######################################')
