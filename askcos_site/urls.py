from django.conf.urls import patterns, include, url
from django.conf import settings
from django.conf.urls.static import static

from django.contrib import admin
admin.autodiscover()

import main.views

# Static (not good for deployment)
urlpatterns = static(settings.STATIC_URL, document_root = settings.STATIC_ROOT)
urlpatterns += static(settings.MEDIA_URL,  document_root = settings.MEDIA_ROOT)

# The rest
urlpatterns += [
    # Examples:
    # url(r'^$', 'askcos_site.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),

    # Admin page
    url(r'^admin/', include(admin.site.urls)),

    # Homepage
    url(r'^$', main.views.index),

    # Retrosynthesis
    url(r'^retro/$', main.views.retro, name = 'retro_home'),
    url(ur'^retro/target=(?P<smiles>.+)$', main.views.retro_target, name = 'retro_target'),

    # Drawing
    url(ur'^draw/smiles/(?P<smiles>.+)$', main.views.draw_smiles, name = 'draw_smiles'),
]