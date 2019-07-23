from django.conf.urls import patterns, include, url
from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
admin.autodiscover()
import django.contrib.auth.urls
from django.views.generic import TemplateView
import askcos_site.main.views as views 
from askcos_site import api

# Static (not good for deployment)
urlpatterns = static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
urlpatterns += static(settings.MEDIA_URL,  document_root=settings.MEDIA_ROOT)

# The rest
urlpatterns += [
    # Examples:
    # url(r'^$', 'askcos_site.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),

    url(r'^robots.txt$', TemplateView.as_view(template_name="robots.txt", content_type="text/plain"), name="robots_file"),
    url(r'^examples/$', TemplateView.as_view(template_name='examples.html')),

    # Admin page
    url(r'^admin/', include(admin.site.urls)),

    # User pages
    url('^registration/', include('registration.backends.simple.urls')),
    # url(r'^login$', views.login, name='login'),
    # url(r'^logout$', views.logout, name='logout'),

    # Homepage
    url(r'^$', views.index, name='index'),
    url(r'^help/modules$', views.modules, name='modules'),
    url(r'^help/tutorial$', views.tutorial, name='tutorial'),
    url(r'^help/faq$', views.faq, name='faq'),

    # Retrosynthesis
    url(r'^retro/$', views.retro, name='retro_home'),
    url(r'^retro/target=(?P<smiles>.+)$', views.retro_target, name='retro_target'),

    # Interactive retrosynthesis
    url(r'^retro_interactive/$', views.retro_interactive, name='retro_interactive'),
    url(r'^retro_interactive/target=(?P<target>.+)$', views.retro_interactive, name='retro_interactive_target'),
    url(r'^retro_interactive_mcts/$', views.retro_interactive_mcts, name='retro_interactive_mcts'),
    url(r'^retro_interactive_mcts/target=(?P<target>.+)$', views.retro_interactive_mcts, name='retro_interactive_mcts_target'),
    url(r'^ajax/smiles_to_image/$', views.ajax_smiles_to_image, name='ajax_smiles_to_image'),
    url(r'^ajax/rxn_to_image/$', views.ajax_rxn_to_image, name='ajax_rxn_to_image'),
    url(r'^ajax/start_retro_celery/$', views.ajax_start_retro_celery, name='ajax_start_retro_celery'),
    url(r'^ajax/start_retro_mcts_celery/$', views.ajax_start_retro_mcts_celery, name='ajax_start_retro_mcts_celery'),
    url(r'^retro_interactive/export/(?P<_id>.+)$', views.export_retro_results, name='export_retro_results'),
    
    # Evaluation
    url(r'^evaluate/$', views.evaluate_rxnsmiles, name='evaluate_rxnsmiles'),
    url(r'^ajax/evaluate_rxnsmiles/$', views.ajax_evaluate_rxnsmiles, name='ajax_evaluate_rxnsmiles'),

    # Context
    url(r'^context/$', views.context_rxnsmiles, name='context_rxnsmiles'),
    url(r'^context/reactants=(?P<reactants>.+)&product=(?P<product>.+)$', views.context_rxnsmiles_target2, name='context_rxnsmiles_target2'),
    url(r'^context/smiles=(?P<smiles>.+)$', views.context_rxnsmiles_target, name='context_rxnsmiles_target'),
    url(r'^ajax/context_rxnsmiles/$', views.ajax_context_rxnsmiles, name='ajax_context_rxnsmiles'),

    # Interactive forward prediction
    url(r'^synth_interactive/$', views.synth_interactive, name='synth_interactive'),
    url(r'^synth_interactive/reactants=(?P<reactants>.+)&reagents=(?P<reagents>.*)&solvent=(?P<solvent>.*)&temperature=(?P<temperature>.*)$', views.synth_interactive, name='synth_interactive_contextspec'),
    url(r'^synth_interactive/reactants=(?P<reactants>.+)&product=(?P<product>.+)$', views.synth_interactive, name='synth_interactive_target2'),
    url(r'^synth_interactive/reactants=(?P<reactants>.+)$', views.synth_interactive, name='synth_interactive_target'),
    url(r'^synth_interactive/smiles=(?P<smiles>.+)$', views.synth_interactive_smiles, name='synth_interactive_target_smiles'),
    url(r'^ajax/start_synth/$', views.ajax_start_synth, name='ajax_start_synth'),
    url(r'^synth_interactive/download$', views.export_synth_results, name='export_synth_results'),

    # Template examination (by str(ObjectID))
    url(r'^template/target=(?P<id>.+)$', views.template_target, name='template_target'),
    url(r'^template/download/target=(?P<id>.+)$', views.template_target_export, name='template_target_export'),
    
    # Reaction examination
    url(r'^reaxys/rxid=(?P<rxid>.+)$', views.rxid_target, name='rxid_target'),

    # Historians
    # url(ur'^history/chemicals/(?P<smiles>.+)$', views.chemical_history_check, name='chemical_history'),
    # url(ur'^history/reactions/(?P<smiles>.+)$', views.reaction_history_check, name='reaction_history'),

    # Pricing
    url(r'^price/$', views.pricing, name='pricing'),
    url(r'^ajax/price_smiles/$', views.ajax_price_smiles, name='ajax_price_smiles'),
    url(r'^price/smiles/(?P<smiles>.+)$', views.price_smiles, name='price_smiles'),
    url(r'^price/xrn/(?P<xrn>.+)$', views.price_xrn, name='price_xrn'),

    # SCScore
    url(r'^scscore/$', views.scscoring, name='scscoring'),
    url(r'^ajax/scscore_smiles/$', views.ajax_scscore_smiles, name='ajax_scscore_smiles'),

    # Drawing
    url(r'^draw/$', views.draw, name='draw'),
    url(r'^draw/smiles/(?P<smiles>.+)$', views.draw_smiles, name='draw_smiles'),
    url(r'^draw/template/(?P<template>.+)$', views.draw_template, name='draw_template'),
    url(r'^draw/reaction/(?P<smiles>.+)$', views.draw_reaction, name='draw_reaction'),
    url(r'^draw/highlight/smiles=(?P<smiles>.+)&reacting_atoms=(?P<reacting_atoms>.+)&bonds=(?P<bonds>.+)$', views.draw_smiles_highlight, name='draw_highlight'),

    # Separation
    # url(r'^separation/input/$', views.sep_input, name='sep_input'),
    url(r'^separation/draw/(?P<fig>.+)$', views.draw_fig, name='draw_fig'),

    # Nearest Neighbor Setup
    # url(r'^nnRecommendation/setup/$', views.nn_predictor_setup, name='setup'),

    # Saved data
    url(r'^saved/$', views.user_saved_results, name='user_saved_results'),
    url(r'^saved/id=(?P<_id>.+)$', views.user_saved_results_id, name='user_saved_results_id'),
    url(r'^ajax/user_save_page/$', views.ajax_user_save_page, name='ajax_user_save_page'),
    url(r'^saved/delete/id=(?P<_id>.+)$', views.user_saved_results_del, name='user_saved_results_del'),

    # Blacklisted reactions
    url(r'^blacklisted/reactions/$', views.user_blacklisted_reactions, name='user_blacklisted_reactions'),
    url(r'^ajax/user_blacklist_reaction/$', views.ajax_user_blacklist_reaction, name='ajax_user_blacklist_reaction'),
    url(r'^blacklisted/reactions/delete/id=(?P<_id>.+)$', views.user_blacklisted_reactions_del, name='user_blacklisted_reactions_del'),
    url(r'^ajax/user_deactivate_reaction/$', views.ajax_user_deactivate_reaction, name='ajax_user_deactivate_reaction'),
    url(r'^ajax/user_activate_reaction/$', views.ajax_user_activate_reaction, name='ajax_user_activate_reaction'),

    # Blacklisted chemicals
    url(r'^blacklisted/chemicals/$', views.user_blacklisted_chemicals, name='user_blacklisted_chemicals'),
    url(r'^ajax/user_blacklist_chemical/$', views.ajax_user_blacklist_chemical, name='ajax_user_blacklist_chemical'),
    url(r'^blacklisted/chemicals/delete/id=(?P<_id>.+)$', views.user_blacklisted_chemicals_del, name='user_blacklisted_chemicals_del'),
    url(r'^ajax/user_deactivate_chemical/$', views.ajax_user_deactivate_chemical, name='ajax_user_deactivate_chemical'),
    url(r'^ajax/user_activate_chemical/$', views.ajax_user_activate_chemical, name='ajax_user_activate_chemical'),
    
    # API endpoints
    url(r'^api/retro/$', api.retro.singlestep, name='retro_api'),
    url(r'^api/fast-filter/$', api.fast_filter.fast_filter, name='fast_filter_api'),
    url(r'^api/context/$', api.context.neural_network, name='context_api'),
    url(r'^api/forward/$', api.forward.template_free, name='forward_api'),
    url(r'^api/template/$', api.template.template, name='template_api'),
    url(r'^api/treebuilder/$', api.tree_builder.tree_builder, name='tree_builder_api'),
    url(r'^api/scscore/$', api.scscore.scscore, name='scscore_api'),
    url(r'^api/price/$', api.price.price, name='price_api'),
    url(r'^api/celery/$', api.celery.celery_status, name='celery_api'),
  
    # Reaction network
    url(r'^retro/network/$', views.retro_network, name='retro_network'),
    
    # Celery status
    url(r'^status/$', views.status),
]
