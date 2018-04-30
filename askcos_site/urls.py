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
    url(r'^$', main.views.index, name='index'),
    url(ur'^modules/$', main.views.modules, name='modules'),

    # Retrosynthesis
    url(r'^retro/$', main.views.retro, name='retro_home'),
    url(ur'^retro/target=(?P<smiles>.+)$', main.views.retro_target, name='retro_target'),

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

    # Context
    url(ur'^context/$', main.views.context_rxnsmiles, name='context_rxnsmiles'),
    url(ur'^context/reactants=(?P<reactants>.+)&product=(?P<product>.+)$', main.views.context_rxnsmiles_target2, name='context_rxnsmiles_target2'),
    url(ur'^context/smiles=(?P<smiles>.+)$', main.views.context_rxnsmiles_target, name='context_rxnsmiles_target'),
    url(ur'^ajax/context_rxnsmiles/$', main.views.ajax_context_rxnsmiles, name='ajax_context_rxnsmiles'),

    # Interactive forward prediction
    url(ur'^synth_interactive/$', main.views.synth_interactive, name='synth_interactive'),
    url(ur'^synth_interactive/reactants=(?P<reactants>.+)&reagents=(?P<reagents>.*)&solvent=(?P<solvent>.*)&temperature=(?P<temperature>.*)$', main.views.synth_interactive, name='synth_interactive_contextspec'),
    url(ur'^synth_interactive/reactants=(?P<reactants>.+)&product=(?P<product>.+)$', main.views.synth_interactive, name='synth_interactive_target2'),
    url(ur'^synth_interactive/reactants=(?P<reactants>.+)$', main.views.synth_interactive, name='synth_interactive_target'),
    url(ur'^synth_interactive/smiles=(?P<smiles>.+)$', main.views.synth_interactive_smiles, name='synth_interactive_target_smiles'),
    url(ur'^ajax/start_synth/$', main.views.ajax_start_synth, name='ajax_start_synth'),
    url(ur'^synth_interactive/download$', main.views.export_synth_results, name='export_synth_results'),

    # Template examination (by str(ObjectID))
    url(r'^template/target=(?P<id>.+)$', main.views.template_target, name='template_target'),

    # Reaction examination
    url(r'^reaxys/rxid=(?P<rxid>.+)$', main.views.rxid_target, name='rxid_target'),

    # Historians
    # url(ur'^history/chemicals/(?P<smiles>.+)$', main.views.chemical_history_check, name='chemical_history'),
    # url(ur'^history/reactions/(?P<smiles>.+)$', main.views.reaction_history_check, name='reaction_history'),

    # Pricing
    url(ur'^price/$', main.views.pricing, name='pricing'),
    url(ur'^ajax/price_smiles/$', main.views.ajax_price_smiles, name='ajax_price_smiles'),
    url(ur'^price/smiles/(?P<smiles>.+)$', main.views.price_smiles, name='price_smiles'),
    url(ur'^price/xrn/(?P<xrn>.+)$', main.views.price_xrn, name='price_xrn'),

    # Drawing
    url(ur'^draw/$', main.views.draw, name='draw'),
    url(ur'^draw/smiles/(?P<smiles>.+)$', main.views.draw_smiles, name='draw_smiles'),
    url(ur'^draw/template/(?P<template>.+)$', main.views.draw_template, name='draw_template'),
    url(ur'^draw/reaction/(?P<smiles>.+)$', main.views.draw_reaction, name='draw_reaction'),

    # Separation
    # url(r'^separation/input/$', main.views.sep_input, name='sep_input'),
    url(r'^separation/draw/(?P<fig>.+)$', main.views.draw_fig, name='draw_fig'),

    # Nearest Neighbor Setup
    # url(r'^nnRecommendation/setup/$', main.views.nn_predictor_setup, name='setup'),

    # Saved data
    url(r'^saved/$', main.views.user_saved_results, name='user_saved_results'),
    url(r'^saved/id=(?P<_id>.+)$', main.views.user_saved_results_id, name='user_saved_results_id'),
    url(ur'^ajax/user_save_page/$', main.views.ajax_user_save_page, name='ajax_user_save_page'),
    url(ur'^saved/delete/id=(?P<_id>.+)$', main.views.user_saved_results_del, name='user_saved_results_del'),

    # Blacklisted reactions
    url(r'^blacklisted/reactions/$', main.views.user_blacklisted_reactions, name='user_blacklisted_reactions'),
    url(ur'^ajax/user_blacklist_reaction/$', main.views.ajax_user_blacklist_reaction, name='ajax_user_blacklist_reaction'),
    url(ur'^blacklisted/reactions/delete/id=(?P<_id>.+)$', main.views.user_blacklisted_reactions_del, name='user_blacklisted_reactions_del'),
    url(ur'^ajax/user_deactivate_reaction/$', main.views.ajax_user_deactivate_reaction, name='ajax_user_deactivate_reaction'),
    url(ur'^ajax/user_activate_reaction/$', main.views.ajax_user_activate_reaction, name='ajax_user_activate_reaction'),

    # Blacklisted chemicals
    url(r'^blacklisted/chemicals/$', main.views.user_blacklisted_chemicals, name='user_blacklisted_chemicals'),
    url(ur'^ajax/user_blacklist_chemical/$', main.views.ajax_user_blacklist_chemical, name='ajax_user_blacklist_chemical'),
    url(ur'^blacklisted/chemicals/delete/id=(?P<_id>.+)$', main.views.user_blacklisted_chemicals_del, name='user_blacklisted_chemicals_del'),
    url(ur'^ajax/user_deactivate_chemical/$', main.views.ajax_user_deactivate_chemical, name='ajax_user_deactivate_chemical'),
    url(ur'^ajax/user_activate_chemical/$', main.views.ajax_user_activate_chemical, name='ajax_user_activate_chemical'),
]