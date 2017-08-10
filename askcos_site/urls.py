from django.conf.urls import patterns, include, url
from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
admin.autodiscover()
import django.contrib.auth.urls

import main.views

# Static (not good for deployment)
urlpatterns = static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
urlpatterns += static(settings.MEDIA_URL,  document_root=settings.MEDIA_ROOT)

# The rest
urlpatterns += [
    # Examples:
    # url(r'^$', 'askcos_site.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),

    # Admin page
    url(r'^admin/', include(admin.site.urls)),

    # User pages
    url('^registration/', include('registration.urls')),
    url(r'^login$', main.views.login, name='login'),
    url(r'^logout$', main.views.logout, name='logout'),

    # Homepage
    url(r'^$', main.views.index),

    # Retrosynthesis
    url(r'^retro/$', main.views.retro, name='retro_home'),
    url(ur'^retro/target=(?P<smiles>.+)$', main.views.retro_target, name='retro_target'),
    url(ur'^retro/lit&target=(?P<smiles>.+)$', main.views.retro_lit_target, name='retro_lit_target'),

    # Interactive retrosynthesis
    url(ur'^retro_interactive/$', main.views.retro_interactive, name='retro_interactive'),
    url(ur'^retro_interactive/target=(?P<target>.+)$', main.views.retro_interactive, name='retro_interactive_target'),
    url(ur'^ajax/smiles_to_image/$', main.views.ajax_smiles_to_image, name='ajax_smiles_to_image'),
    url(ur'^ajax/rxn_to_image/$', main.views.ajax_rxn_to_image, name='ajax_rxn_to_image'),
    url(ur'^ajax/start_retro_celery/$', main.views.ajax_start_retro_celery, name='ajax_start_retro_celery'),
    url(ur'^retro_interactive/export/(?P<_id>.+)$', main.views.export_retro_results, name='export_retro_results'),
    
    # Evaluation
    url(ur'^evaluate/$', main.views.evaluate_rxnsmiles, name='evaluate_rxnsmiles'),
    url(ur'^ajax/evaluate_rxnsmiles/$', main.views.ajax_evaluate_rxnsmiles, name='ajax_evaluate_rxnsmiles'),

    # Interactive forward prediction
    url(ur'^synth_interactive/$', main.views.synth_interactive, name='synth_interactive'),
    url(ur'^synth_interactive/reactants=(?P<reactants>.+)&product=(?P<product>.+)$', main.views.synth_interactive, name='synth_interactive_target2'),
    url(ur'^synth_interactive/reactants=(?P<reactants>.+)$', main.views.synth_interactive, name='synth_interactive_target'),
    url(ur'^synth_interactive/smiles=(?P<smiles>.+)$', main.views.synth_interactive_smiles, name='synth_interactive_target_smiles'),
    url(ur'^ajax/start_synth/$', main.views.ajax_start_synth, name='ajax_start_synth'),
    url(ur'^synth_interactive/download$', main.views.export_synth_results, name='export_synth_results'),

    # Forward synthesis
    # url(r'^synth/$', main.views.synth, name='synth_home'),
    # url(ur'^synth/target=(?P<smiles>.+)$', main.views.synth_target, name='synth_target'),

    # Template examination (by str(ObjectID))
    url(r'^template/target=(?P<id>.+)$', main.views.template_target, name='template_target'),

    # Reaction examination
    url(r'^reaxys/rxid=(?P<rxid>.+)$', main.views.rxid_target, name='rxid_target'),

    # Pricing
    url(ur'^price/$', main.views.pricing, name='pricing'),
    url(ur'^ajax/price_smiles/$', main.views.ajax_price_smiles, name='ajax_price_smiles'),
    url(ur'^price/smiles/(?P<smiles>.+)$', main.views.price_smiles, name='price_smiles'),
    url(ur'^price/xrn/(?P<xrn>.+)$', main.views.price_xrn, name='price_xrn'),

    # Drawing
    url(ur'^draw/$', main.views.draw, name='draw'),
    url(ur'^draw/smiles/(?P<smiles>.+)$', main.views.draw_smiles, name='draw_smiles'),
    url(ur'^draw/smiles_page/(?P<smiles>.+)$', main.views.draw_smiles_page, name='draw_smiles_page'),
    url(ur'^draw/template/(?P<template>.+)$', main.views.draw_template, name='draw_template'),
    url(ur'^draw/template_page/(?P<template>.+)$', main.views.draw_template_page, name='draw_template_page'),
    url(ur'^draw/reaction/(?P<smiles>.+)$', main.views.draw_reaction, name='draw_reaction'),
    url(ur'^draw/reaction_page/(?P<smiles>.+)$', main.views.draw_reaction_page, name='draw_reaction_page'),

    # Showing a path (testing)
    url(ur'^draw/synthesis_tree/$', main.views.draw_synthesis_tree, name='draw_synthesis_tree'),
    url(ur'^draw/synthesis_tree/id=(?P<id>.+)$', main.views.draw_synthesis_tree, name='draw_synthesis_tree_click'),

    # Separation
    url(r'^separation/input/$', main.views.sep_input, name='sep_input'),
    url(r'^separation/draw/(?P<fig>.+)$', main.views.draw_fig, name='draw_fig'),

    # Nearest Neighbor Setup
    url(r'^nnRecommendation/setup/$', main.views.nn_predictor_setup, name='setup'),

    # Saved data
    url(r'^saved/$', main.views.user_saved_results, name='user_saved_results'),
    url(r'^saved/id=(?P<_id>.+)$', main.views.user_saved_results_id, name='user_saved_results_id'),
    url(ur'^ajax/user_save_page/$', main.views.ajax_user_save_page, name='ajax_user_save_page'),
    url(ur'^saved/delete/id=(?P<_id>.+)$', main.views.user_saved_results_del, name='user_saved_results_del'),

    # Blacklisted reactions
    url(r'^blacklisted/$', main.views.user_blacklisted_reactions, name='user_blacklisted_reactions'),
    url(ur'^ajax/user_blacklist_reaction/$', main.views.ajax_user_blacklist_reaction, name='ajax_user_blacklist_reaction'),
    url(ur'^blacklisted/delete/id=(?P<_id>.+)$', main.views.user_blacklisted_reactions_del, name='user_blacklisted_reactions_del'),
    url(ur'^ajax/user_deactivate_reaction/$', main.views.ajax_user_deactivate_reaction, name='ajax_user_deactivate_reaction'),
    url(ur'^ajax/user_activate_reaction/$', main.views.ajax_user_activate_reaction, name='ajax_user_activate_reaction'),
]