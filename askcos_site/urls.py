from django.conf.urls import patterns, include, url
from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
admin.autodiscover()
import django.contrib.auth.urls

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

    # User pages
    url(r'^login$', main.views.login, name = 'login'),
    url(r'^logout$', main.views.logout, name = 'logout'),

    # Homepage
    url(r'^$', main.views.index),

    # Retrosynthesis
    url(r'^retro/$', main.views.retro, name = 'retro_home'),
    url(ur'^retro/target=(?P<smiles>.+)$', main.views.retro_target, name = 'retro_target'),

    # Forward synthesis
    url(r'^synth/$', main.views.synth, name = 'synth_home'),
    url(ur'^synth/target=(?P<smiles>.+)$', main.views.synth_target, name = 'synth_target'),

    # Template examination (by str(ObjectID))
    url(r'^template/target=(?P<id>.+)$', main.views.template_target, name = 'template_target'),

    # Drawing
    url(ur'^draw/$', main.views.draw, name = 'draw'),
    url(ur'^draw/smiles/(?P<smiles>.+)$', main.views.draw_smiles, name = 'draw_smiles'),
    url(ur'^draw/smiles_page/(?P<smiles>.+)$', main.views.draw_smiles_page, name = 'draw_smiles_page'),
    url(ur'^draw/template/(?P<template>.+)$', main.views.draw_template, name = 'draw_template'),
    url(ur'^draw/template_page/(?P<template>.+)$', main.views.draw_template_page, name = 'draw_template_page'),
    url(ur'^draw/reaction/(?P<smiles>.+)$', main.views.draw_reaction, name = 'draw_reaction'),
    url(ur'^draw/reaction_page/(?P<smiles>.+)$', main.views.draw_reaction_page, name = 'draw_reaction_page'),
]