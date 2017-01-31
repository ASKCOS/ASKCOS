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
    url(ur'^retro/lit&target=(?P<smiles>.+)$', main.views.retro_lit_target, name = 'retro_lit_target'),

    # Interactive retrosynthesis
    url(ur'^retro_interactive/$', main.views.retro_interactive, name = 'retro_interactive'),
    url(ur'^ajax/smiles_to_image_retro/$', main.views.ajax_smiles_to_image_retro, name = 'ajax_smiles_to_image_retro'),
    url(ur'^ajax/start_retro/$', main.views.ajax_start_retro, name = 'ajax_start_retro'),
    url(ur'^ajax/pause_retro/$', main.views.ajax_pause_retro, name = 'ajax_pause_retro'),
    url(ur'^ajax/stop_retro/$', main.views.ajax_stop_retro, name = 'ajax_stop_retro'),
    url(ur'^ajax/update_retro_stats/$', main.views.ajax_update_retro_stats, name = 'ajax_update_retro_stats'),
    url(ur'^ajax/update_retro/$', main.views.ajax_update_retro, name = 'ajax_update_retro'),

    # Interactive forward prediction
    url(ur'^synth_interactive/$', main.views.synth_interactive, name = 'synth_interactive'),
    url(ur'^ajax/smiles_to_image_synth/$', main.views.ajax_smiles_to_image_synth, name = 'ajax_smiles_to_image_synth'),
    url(ur'^ajax/start_synth/$', main.views.ajax_start_synth, name = 'ajax_start_synth'),

    # Forward synthesis
    url(r'^synth/$', main.views.synth, name = 'synth_home'),
    url(ur'^synth/target=(?P<smiles>.+)$', main.views.synth_target, name = 'synth_target'),

    # Template examination (by str(ObjectID))
    url(r'^template/target=(?P<id>.+)$', main.views.template_target, name = 'template_target'),

    # Reaction examination
    url(r'^reaxys/rxid=(?P<rxid>.+)$', main.views.rxid_target, name = 'rxid_target'),

    # Pricing
    url(ur'^price/smiles/(?P<smiles>.+)$', main.views.price_smiles, name = 'price_smiles'),
    url(ur'^price/xrn/(?P<xrn>.+)$', main.views.price_xrn, name = 'price_xrn'),

    # Drawing
    url(ur'^draw/$', main.views.draw, name = 'draw'),
    url(ur'^draw/smiles/(?P<smiles>.+)$', main.views.draw_smiles, name = 'draw_smiles'),
    url(ur'^draw/smiles_page/(?P<smiles>.+)$', main.views.draw_smiles_page, name = 'draw_smiles_page'),
    url(ur'^draw/template/(?P<template>.+)$', main.views.draw_template, name = 'draw_template'),
    url(ur'^draw/template_page/(?P<template>.+)$', main.views.draw_template_page, name = 'draw_template_page'),
    url(ur'^draw/reaction/(?P<smiles>.+)$', main.views.draw_reaction, name = 'draw_reaction'),
    url(ur'^draw/reaction_page/(?P<smiles>.+)$', main.views.draw_reaction_page, name = 'draw_reaction_page'),

    # Showing a path (testing)
    url(ur'^draw/synthesis_tree/$', main.views.draw_synthesis_tree, name = 'draw_synthesis_tree'),
    url(ur'^draw/synthesis_tree/id=(?P<id>.+)$', main.views.draw_synthesis_tree, name = 'draw_synthesis_tree_click'),

    # Using Ajax calls (testing)
    url(ur'^test/the_time/$', main.views.the_time, name = 'the_time'),
    url(ur'^ajax/get_the_time/$', main.views.get_the_time, name = 'get_the_time'),
]