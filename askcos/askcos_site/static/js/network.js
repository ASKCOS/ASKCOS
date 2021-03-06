var container = document.getElementsByClassName('container')[0];
container.classList.remove('container')
container.classList.add('container-fluid')
container.style.width=null;

function num2str(n, len) {
    if (len == undefined) {
        return n == undefined || isNaN(n) ? 'N/A' : n.toString();
    } else {
        return n == undefined || isNaN(n) ? 'N/A' : n.toFixed(len);
    }
}

// check whether the set on which the  
// method is invoked is the subset of  
// otherset or not 
function subSet(s, otherSet) {
    if(s.size > otherSet.size) {
        return false;
    } else {
        for(var elem of s) {
            // if any of the element of
            // this is not present in the
            // otherset then return false
            if(!otherSet.has(elem))
                return false;
        }
        return true;
    }
};

function getCookie(cname) {
    var name = cname + "=";
    var cookie_str = document.cookie;
    if (cookie_str && cookie_str != '') {
        var cookie_splitted = cookie_str.split(';');
        for(var i = 0; i <cookie_splitted.length; i++) {
            var c = cookie_splitted[i].trim();
            if (c.indexOf(name) == 0) {
                return decodeURIComponent(c.substring(name.length, c.length));
            }
        }
    }
  return undefined;
}

function copyToClipboard(text) {
    var dummy = document.createElement("textarea");
    document.body.appendChild(dummy);
    dummy.value = text;
    dummy.select();
    document.execCommand("copy");
    document.body.removeChild(dummy);
}

function showLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "block";
}

function hideLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "none";
}

function ctaToNode(cta, id) {
    if (cta.is_reaction) {
        return {
            id: id,
            graphId: cta.id,
            plausibility: cta.plausibility,
            numExamples: cta.num_examples,
            reactionSmiles: cta.smiles,
            templateIds: cta.tforms,
            templateScore: cta.template_score,
            retroscore: cta.template_score,
            type: 'reaction',
        }
    }
    else {
        let node = {
            id: id,
            graphId: cta.id,
            asProduct: cta.as_product,
            asReactant: cta.as_reactant,
            type: 'chemical',
            ppg: cta.ppg,
            smiles: cta.smiles,
            image: "/draw/smiles/"+encodeURIComponent(cta.smiles),
            shape: 'image',
            borderWidth: 2
        }
        if (node.ppg == 0) {
            node.color = {border: '#880000'}
        }
        else {
            node.color = {border: '#008800'}
        }
        return node
    }
}

function addReactions(reactions, sourceNode, nodes, edges, reactionLimit) {
    var added = 0
    for (r of reactions) {
        if (added >= reactionLimit) {
            break;
        }
        if (r.show) {
            addReaction(r, sourceNode, nodes, edges);
            added += 1;
        }
    }
}

function addReaction(reaction, sourceNode, nodes, edges) {
    var rId = nodes.max('id').id+1;
    nodes.add({
        id: rId,
        label: '#'+reaction['rank'],
        rank: reaction['rank'],
        ffScore: num2str(reaction['plausibility'] ,3),
        retroscore: num2str(reaction['score'], 3),
        templateScore: num2str(reaction['template_score'], 3),
        numExamples: num2str(reaction['num_examples']),
        templateIds: reaction['templates'],
        reactionSmiles: reaction.smiles+'>>'+sourceNode.smiles,
        type: 'reaction',
        value: 1,
        mass: 1
    })
    if (edges.max('id')) {
        var eId = edges.max('id').id+1
    }
    else {
        var eId = 0
    }
    edges.add({
        id: eId,
        from: sourceNode.id,
        to: rId,
                scaling: {
                    min: 1,
                    max: 5,
                    customScalingFunction: function(min, max, total, value) {
                        if (value > 0.25) {
                            return 1.0
                        }
                        else{
                            return 16*value*value
                        }
                    }
                },
        color: {
            color: '#000000',
            inherit: false
        },
        value: Number(reaction['template_score'])
    })
    for (n in reaction['smiles_split']) {
        var smi = reaction['smiles_split'][n];
        fetch('/api/buyables/search/?q='+encodeURIComponent(smi)+'&canonicalize=True')
        .then(resp => resp.json())
        .then(json => {
            var mysmi = json['search'];
            if (json.buyables.length > 0) {
                var ppg = json.buyables[0].ppg
                var buyable = true
                var source = json.buyables[0].source
            }
            else {
                var ppg = "not buyable"
                var buyable = false
                var source = ''
            }
            if (buyable) {
                var color = "#008800"
            }
            else {
                var color = "#880000"
            }
            var nId = nodes.max('id').id+1;
            nodes.add({
                id: nId,
                smiles: mysmi,
                image: app.getMolDrawEndPoint(mysmi),
                shape: "image",
                borderWidth: 2,
                type: 'chemical',
                mass: 1,
                value: 10,
                ppg: ppg,
                source: source,
                color: {
                    border: color
                }
            })
            edges.add({
                id: edges.max('id').id+1,
                from: rId,
                to: nId,
                scaling: {
                    min: 1,
                    max: 5,
                    customScalingFunction: function(min, max, total, value) {
                        if (value > 0.25) {
                            return 1.0
                        }
                        else{
                            return 16*value*value
                        }
                    }
                },
                color: {
                    color: '#000000',
                    inherit: false
                },
                value: Number(reaction['template_score'])
            })
        })
        reaction.inViz = true;
    }
}

function parentOf(id, nodes, edges) {
    var parentId = -1;
    edges.forEach(function(e) {
        if (e!=null && e.to==id) {
            parentId = e.from
        }
    })
    return parentId
}

function childrenOf(id, nodes, edges) {
    var children = [];
    edges.forEach(function(e) {
        if (e!=null && e.from==id) {
            children.push(e.to)
        }
    })
    return children
}

function allChildrenOf(id, nodes, edges) {
    var children = [];
    edges.forEach(function(e) {
        if (e!=null && e.from==id) {
            children.push(e.to);
            var tmpChildren = allChildrenOf(e.to, nodes, edges);
            for (n in tmpChildren) {
                var child = tmpChildren[n];
                children.push(child);
            }
        }
    })
    return children
}

function removeChildrenFrom(id, nodes, edges) {
    var children = allChildrenOf(id, nodes, edges);
    nodes.remove(children);
}

function cleanUpEdges(nodes, edges) {
    var nodeIds = nodes.getIds();
    edges.forEach(function(edge) {
        if (!nodeIds.includes(edge.from) | !nodeIds.includes(edge.to)) {
            edges.remove(edge.id)
        }
    })
}

var network;

function initializeNetwork(data, hierarchical) {
    var container = document.getElementById('network');
    var options = {
        nodes: {
            color: {
                border: '#000000',
                background: '#FFFFFF'
            },
            shapeProperties: {
                useBorderWithImage:true
            }
        },
        edges: {
            length: 1
        },
        interaction: {
            multiselect: true
        },
    };
    if (hierarchical) {
        options.layout = {
          hierarchical: {
            levelSeparation: 150,
            nodeSpacing: 175,
          }
        }
        options.physics = {
            stabilization: false,
            barnesHut: {
              gravitationalConstant: -80000,
            //   springConstant: 0.001,
            //   springLength: 200
            }
        }
    }
    network = new vis.Network(container, data, options);
    network.on("beforeDrawing",  function(ctx) {
        ctx.save();
        ctx.setTransform(1, 0, 0, 1, 0, 0);
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height)
        ctx.restore();
    })
    return network
}

/* DnD */
function disable_dragstart_handler(e) {
    e.preventDefault();
}

function clusteredit_dragover_handler(event) {
    event.preventDefault(); // important
    event.stopPropagation();
    event.dataTransfer.dropEffect = 'move'; // important
}

function clusteredit_dragenter_handler(event) {
    event.preventDefault(); // important
    event.stopPropagation();
    event.target.classList.add('dragover');
}

function clusteredit_dragleave_handler(event) {
    event.target.classList.remove('dragover');
}

Vue.component('modal', {
    template: '#modal-template'
})

var app = new Vue({
    el: '#app',
    data: {
        window: {
            width: 0,
            height: 0,
        },
        target: '',
        network: {},
        data: {
            nodes: {},
            edges: {}
        },
        results: {},
        templateNumExamples: {},
        nodeStructure: {},
        allowCluster: true,
        allowResolve: false,
        showSettingsModal: false,
        showLoadModal: false,
        showDownloadModal: false,
        showClusterPopoutModal: false,
        showClusterEditModal: false,
        showAddNewPrecursorModal: false,
        downloadName: "network.json",
        modalData: {},
        clusterPopoutModalData: {
            optionsDisplay : {
                showScore: false,
                showNumExample: true,
                showTemplateScore: false,
                showPlausibility: true,
                showClusterId: false,
            },
        },
        clusterEditModalData: {
            optionsDisplay : {
                showScore: false,
                showNumExample: false,
                showTemplateScore: false,
                showPlausibility: false,
                showClusterId: false,
            },
        },
        addNewPrecursorModal: {},
        clusterOptions: {
            allowRemovePrecursor: true,
            feature: 'original',
            fingerprint:'morgan',
            fpRadius: 1, fpBits: 512,
            cluster_method: 'kmeans',
            isAlternatingColor: false,
        },
        selected: null,
        isHighlightAtom: true,
        reactionLimit: 5,
        templatePrioritization: "reaxys",
        templateSet: "reaxys",
        precursorScoring: "RelevanceHeuristic",
        numTemplates: 1000,
        maxCumProb: 0.999,
        minPlausibility: 0.01,
        sortingCategory: "score",
        networkHierarchical: false
    },
    beforeMount: function() {
        this.allowResolve = this.$el.querySelector('[ref="allowResolve"]').checked;
    },
    created: function() {
        window.addEventListener('resize', this.handleResize);
        this.handleResize();
    },
    mounted: function() {
        var urlParams = new URLSearchParams(window.location.search);
        let loadTreeBuilder = urlParams.get('tb')
        let numTrees = urlParams.get('view')
        console.log(loadTreeBuilder)
        if (loadTreeBuilder) {
            this.loadFromTreeBuilder(loadTreeBuilder, numTrees)
        }
    },
    destroyed: function() {
        window.removeEventListener('resize', this.handleResize);
    },
    methods: {
        handleResize: function() {
            this.window.width = window.innerWidth;
            this.window.height = window.innerHeight;
        },
        requestUrl: function(smiles) {
            var url = '/api/retro/?';
            var params = {
                target: smiles,
                template_set: this.templateSet,
                template_prioritizer: this.templatePrioritization,
                precursor_prioritization: this.precursorScoring,
                num_templates: this.numTemplates,
                max_cum_prob: this.maxCumProb,
                filter_threshold: this.minPlausibility,
                cluster_method: this.clusterOptions.cluster_method,
                cluster_feature: this.clusterOptions.feature,
                cluster_fp_type: this.clusterOptions.fingerprint,
                cluster_fp_length: this.clusterOptions.fpBits,
                cluster_fp_radius: this.clusterOptions.fpRadius
            }
            var queryString = Object.keys(params).map((key) => {
                return encodeURIComponent(key) + '=' + encodeURIComponent(params[key])
            }).join('&');
            return url+queryString;
        },
        resolveChemName: function(name) {
            if (this.allowResolve) {
                var url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+encodeURIComponent(name)+'/property/IsomericSMILES/txt'
                console.log(url)
                var text = fetch(url)
                    .then(resp => {
                        if (resp.status == 404) {
                            throw Error(resp.statusText);
                        } else {
                            return resp.text()
                        }
                    })
                    .catch(err => {
                        throw Error('Cannot resolve "'+name+'" to smiles: '+err.message);
                    })
                return text;
            } else {
                throw Error('Resolving chemical name using external server is not allowed.');
            }
        },
        validatesmiles: function(s, iswarning) {
            var url = '/api/validate-chem-name/?smiles='+encodeURIComponent(s)
            console.log(url)
            var res = fetch(url)
                .then(resp => {
                    if (!resp.ok) {
                        throw 'Unable to connect to server: error code '+resp.status;
                    } else {
                        return resp.json()
                    }
                })
                .then(res_json => {
                    if (!res_json['correct_syntax']) {
                        if(iswarning) alert('Input SMILES string: invalid syntax')
                        return false
                    } else if (!res_json['valid_chem_name']) {
                        if(iswarning) alert('Input SMILES string: invalid chemical')
                        return false
                    } else {
                        return true
                    }
                })
            return res;
        },
        changeTarget: function() {
            showLoader();
            this.validatesmiles(this.target, !this.allowResolve)
            .then(isvalidsmiles => {
                if (isvalidsmiles) {
                    return this.target
                } else {
                    return this.resolveChemName(this.target)
                }
            })
            .then(x => {
                this.target = x;
                if (this.target != undefined) {
                    var url = this.requestUrl(this.target);
                    fetch(url)
                    .then(resp => {
                        if (!resp.ok) {
                            throw resp.status;
                        }
                        return resp;
                    })
                    .then(resp => resp.json())
                    .then(json => {
                        if ('error' in json) {
                            throw json['error']
                        } else {
                            this.data.nodes = new vis.DataSet([
                                this.createTargetNode(this.target)
                            ])
                            this.data.edges = new vis.DataSet([]);
                            initializeNetwork(this.data);
                            network.on('selectNode', this.showInfo);
                            network.on('deselectNode', this.clearSelection);
                            this.$set(this.results, this.target, json['precursors']);
                            this.initClusterShowCard(this.target); // must be called immediately after adding results
                            addReactions(this.results[this.target], this.data.nodes.get(0), this.data.nodes, this.data.edges, this.reactionLimit);
                            this.getTemplateNumExamples(this.results[this.target]);
                            hideLoader();
                            fetch('/api/buyables/search/?q='+encodeURIComponent(this.target)+'&canonicalize=True')
                                .then(resp => resp.json())
                                .then(json => {
                                    if (json.buyables.length > 0) {
                                        var ppg = json.buyables[0].ppg
                                    }
                                    else {
                                        var ppg = "not buyable"
                                    }
                                    this.data.nodes.update({id: 0, ppg: ppg});
                                    network.selectNodes([0]);
                                    this.selected = this.data.nodes.get(0);
                            })
                        }
                    })
                    .catch(error => {
                        hideLoader();
                        alert('There was an error fetching precursors for this target with the supplied settings: '+error)
                    })
                } else {
                    hideLoader();
                }
            })
            .catch(error => {
                hideLoader();
                var error_msg = 'unknown error'
                if ('message' in error) {
                    error_msg = error.name+':'+error.message
                } else if (typeof(error) == 'string') {
                    error_msg = error
                }
                alert('There was an error fetching precursors for this target with the supplied settings: '+error_msg)
            })
        },
        toggleHierarchical: function() {
          if (typeof(network) == 'undefined') {
            return
          }
          if (this.networkHierarchical) {
            network.setOptions({'layout': {'hierarchical': false}})
            document.querySelector('#hierarchical-button').innerHTML = 'G'
            this.networkHierarchical = false;
          }
          else {
            network.setOptions({'layout': {'hierarchical': {sortMethod: 'directed'}}})
            document.querySelector('#hierarchical-button').innerHTML = 'H'
            this.networkHierarchical = true;
          }
        },
        expandNode: function() {
            if (this.isModalOpen() || typeof(network) == "undefined") {
                return
            }
            showLoader();
            var selected = network.getSelectedNodes();
            if (selected.length != 1) {
              hideLoader();
              if (selected.length == 0) {
                  alert('Please select a terminal chemical node to expand')
              }
              else {
                  alert('Please only select 1 node at a time to expand')
              }
              return
            }
            var nodeId = selected[0];
            if (typeof(nodeId) == 'string' && nodeId.startsWith('cluster')) {
                alert('Cannot expand collpased node! To toggle collpased state, click collapse toggle button again with collapsed cluster selected.')
                hideLoader();
                return
            }
            var node = this.data.nodes.get(nodeId)
            if (node.type != 'chemical') {
                alert('Cannot expand reaction; try expanding with a chemical node selected');
                hideLoader();
                return
            }
            var childrenOfSelected = childrenOf(nodeId, this.data.nodes, this.data.edges);
            if (childrenOfSelected.length != 0) {
                alert("You've already expanded this node. If you would like to re-expand, please use the 'Remove children nodes' button to clear results in the visualization for this chemical. Please note that this will replace the previously predicted results for this chemical (for example, if you've changed any settings)")
                hideLoader();
                return
            }
            var smi = node.smiles;
            var url = this.requestUrl(smi);
            fetch(url)
                .then(resp => {
                    if (!resp.ok) {
                        throw Error(resp.statusText);
                    }
                    return resp;
                })
                .then(resp => resp.json())
                .then(json => {
                    var reactions = json['precursors'];
                    this.$set(this.results, smi, reactions);
                    if (reactions.length==0) {
                        alert('No precursors found!')
                    }
                    this.initClusterShowCard(smi); // must be called immediately after adding results
                    addReactions(this.results[smi], this.data.nodes.get(nodeId), this.data.nodes, this.data.edges, this.reactionLimit);
                    this.getTemplateNumExamples(this.results[smi]);
                    this.selected = node;
                    this.reorderResults();
                    hideLoader();
                })
                .catch(error => {
                    hideLoader();
                    alert('There was an error fetching precursors for this target with the supplied settings')
                })
        },
        getTemplateNumExamples: function(reactions) {
            for (reaction of reactions) {
                for (templateId of reaction['templates']) {
                    this.apiTemplateCount(templateId);
                }
            }
        },
        apiTemplateCount: function(templateId) {
            if (typeof(this.templateNumExamples[templateId]) == 'undefined') {
                fetch('/api/template/?id='+templateId)
                .then(resp => resp.json())
                .then(json => {
                    var id = json["request"]["id"][0];
                    var count = json["template"]["count"];
                    this.templateNumExamples[id] = count;
                })
            }
        },
        deleteChoice: function() {
            var selected = network.getSelectedNodes();
            for (n in selected) {
                var nodeId = selected[n];
                if (this.data.nodes.get(nodeId).type=='chemical') {
                    let res = confirm('This will delete all children nodes of the currently selected node. Continue?')
                    if (res) {
                        this.deleteChildren()
                    }
                }
                else {
                    let res = confirm('This will delete the currently selected node and all children node. Continue?')
                    if (res) {
                        this.deleteNode()
                    }
                }
            }
        },
        deleteNode: function() {
            if (this.isModalOpen() || typeof(network) == "undefined") {
                return
            }
            var selected = network.getSelectedNodes();
            for (n in selected) {
                var nodeId = selected[n];
                if (this.data.nodes.get(nodeId).type=='chemical') {
                    alert('We do not allow deleting chemical nodes! It will leave its parent reaction node missing information! Only delete reaction nodes with this button.')
                    continue
                }
                var node = this.data.nodes.get(nodeId);
                var parentNodeId = parentOf(nodeId, this.data.nodes, this.data.edges);
                var parentNode = this.data.nodes.get(parentNodeId);
                for (result of this.results[parentNode.smiles]) {
                    if (result.rank == node.rank) {
                        result.inViz = false;
                        break;
                    }
                }
                removeChildrenFrom(nodeId, this.data.nodes, this.data.edges);
                this.data.nodes.remove(nodeId);
            }
            cleanUpEdges(this.data.nodes, this.data.edges);
        },
        deleteChildren: function() {
            var selected = network.getSelectedNodes();
            for (n in selected) {
                var nodeId = selected[n];
                if (this.data.nodes.get(nodeId).type=='reaction') {
                    alert('We do not allow deleting children of reaction nodes! It will leave the reaction node missing information! Only delete children of chemical nodes with this button.')
                    continue
                }
                var node = this.data.nodes.get(nodeId);
                for (result of this.results[node.smiles]) {
                    result.inViz = false;
                }
                removeChildrenFrom(nodeId, this.data.nodes, this.data.edges);
            }
            cleanUpEdges(this.data.nodes, this.data.edges);
            var oldSelected = this.selected;
            this.selected = null;
            network.unselectAll();
//             this.selected = oldSelected;
        },
        toggleResolver: function() {
            if (this.allowResolve) {
                this.allowResolve = false
            }
            else {
                this.allowResolve = true
            }
        },
        download: function() {
            if (this.data.nodes.length == null) {
                alert("There's no network to download!")
                return
            }
            var downloadData = {nodes: [], edges: [], results: this.results}
            this.data.nodes.forEach(function(e) {
                downloadData.nodes.push(e)
            })
            this.data.edges.forEach(function(e) {
                downloadData.edges.push(e)
            })
            var dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(downloadData));
            var dlAnchorElem = document.getElementById('downloadAnchorElem');
            dlAnchorElem.setAttribute("href",     dataStr     );
            dlAnchorElem.setAttribute("download", this.downloadName);
            dlAnchorElem.click();
        },
        hasUndefinedGroupid: function() {
            // check if this.results has group_id
            for (s in this.results) {
                var precursors = this.results[s];
                for (i of precursors) {
                    if (i.group_id == undefined) {
                        return true;
                    }
                }
            }
            return false;
        },
        load: function() {
            var file = document.getElementById("loadNetwork").files[0];
            var reader = new FileReader();
            var app = this;
            reader.onload = (function(theFile) {return function(e) {
                var data = JSON.parse(e.target.result);
                app.target = data.nodes[0].smiles;
                app.data.nodes = new vis.DataSet(data.nodes);
                app.data.edges = new vis.DataSet(data.edges);
                app.results = data.results;
                if (app.hasUndefinedGroupid()) {
                    let res = confirm('The uploaded json file does not have reaction cluster information for some precursors. Select "OK" to re-cluster all of them. This will erase existing reaction cluster information. Select "Cancel" to skip, however, reaction cluster function may not work correctly until you re-cluster manually.');
                    if (res) {
                        for (s in app.results) {
                            app.requestClusterId(s);
                        }
                    } else {
                        app.allowCluster = false;
                    }
                }
                network = initializeNetwork(app.data)
                network.on('selectNode', app.showInfo);
                network.on('deselectNode', app.clearSelection);
            }})(file);
            reader.readAsText(file)
        },
        clear: function() {
            this.target = '';
            this.data.nodes.remove(this.data.nodes.getIds());
            this.data.edges.remove(this.data.edges.getIds());
            this.selected = null;
            document.querySelector('#hierarchical-button').innerHTML = 'G';
        },
        clearSelection: function() {
            this.selected = null;
        },
        copySelectedSmiles: function() {
            var copyTooltip = document.querySelector('#copy-tooltip')
            if (this.selected.type == 'chemical') {
                copyToClipboard(this.selected.smiles)
            }
            else {
                copyToClipboard(this.selected.reactionSmiles)
            }
            copyTooltip.innerHTML = 'Copied!'
            setTimeout(() => {copyTooltip.innerHTML = "Click to copy!"}, 2000)
        },
        collapseNode: function() {
            var selected = network.getSelectedNodes();
            var type = typeof selected[0];
            if (type == "string") {
                network.openCluster(selected[0])
            }
            else {
                var forCluster = allChildrenOf(selected[0], app.data.nodes, app.data.edges)
                var options = {
                    joinCondition:function(nodeOptions) {
                        if (forCluster.includes(nodeOptions.id) | nodeOptions.id==selected[0]) {
                            return true;
                        }
                    }
                }
                network.clustering.cluster(options);
            }
        },
        addFromResults: function(selected, reaction) {
            if (reaction.inViz) {
                return
            }
            addReaction(reaction, selected, this.data.nodes, this.data.edges);
            reaction.inViz = true;
            document.querySelector('.addRes[data-rank="'+Number(reaction.rank)+'"]').style.display='none';
            document.querySelector('.remRes[data-rank="'+Number(reaction.rank)+'"]').style.display='';
        },
        remFromResults: function(selected, reaction) {
            var rsmi = reaction.smiles+'>>'+selected.smiles;
            var selectedChildren = this.data.nodes.get(childrenOf(selected.id, this.data.nodes, this.data.edges));
            for (var child of selectedChildren) {
                if (child.reactionSmiles == rsmi) {
                    removeChildrenFrom(child.id, this.data.nodes, this.data.edges);
                    this.data.nodes.remove(child.id);
                    cleanUpEdges(this.data.nodes, this.data.edges);
                    document.querySelector('.addRes[data-rank="'+Number(reaction.rank)+'"]').style.display='';
                    document.querySelector('.remRes[data-rank="'+Number(reaction.rank)+'"]').style.display='none';
                    reaction.inViz = false;
                    break;
                }
            }
        },
        resetSortingCategory: function() {
            this.sortingCategory = 'score'
            this.reorderResults()
        },
        reorderResults: function() {
            var sortingCategory = this.sortingCategory;
            if (this.selected.type != 'chemical') {
                return
            }
            var smiles = this.selected.smiles;
            var results = this.results[smiles];
            if (typeof(results) == 'undefined') {
                return
            }
            results.sort((a, b) => {
                var a_ = a[sortingCategory] == undefined ? 0 : a[sortingCategory];
                var b_ = b[sortingCategory] == undefined ? 0 : b[sortingCategory];
                return b_ - a_;
            })
            var prevSelected = this.selected;
            this.selected = undefined;
            this.selected = prevSelected;
        },
        showInfo: function(obj) {
            var nodeId = obj.nodes[obj.nodes.length-1];
            var node = this.data.nodes.get(nodeId);
            if (node == null) {
                return
            }
            this.selected = node;
            this.reorderResults();
            if (node.type == 'chemical') {
                console.log('chemical', node);
                if (typeof(this.results[node.smiles]) != 'undefined') {
                    console.log(this.results[node.smiles])
                }
            }
            else if (node.type == 'reaction') {
                console.log('reaction', node)
            }
        },
        openModal: function(modalName) {
            /*
            this.clearSelection();
            if (network) {
                network.unselectAll();
            }
            */
            if (modalName == "settings") {
                this.showSettingsModal = true
            }
            else if (modalName == "download") {
                this.showDownloadModal = true
            }
            else if (modalName == "load") {
                this.showLoadModal = true
            }
        },
        openClusterPopoutModal: function(selected, res) {
            if(selected == undefined) {
                alert('No target molecule selected. Please select a molecule in the tree.')
                return
            }
            /*
             * cannot deselect
            this.clearSelection();
            if (network) {
                network.unselectAll();
            }
            */
            this.$set(this.clusterPopoutModalData, 'selected', selected);
            this.$set(this.clusterPopoutModalData, 'selectedSmiles', selected.smiles);
            this.$set(this.clusterPopoutModalData, 'res', res);
            this.$set(this.clusterPopoutModalData, 'group_id', res.group_id);
            this.showClusterPopoutModal = true;
        },
        closeClusterPopoutModal: function() {
            this.showClusterPopoutModal = false;
            this.clusterPopoutModalData['selected'] = undefined;
            this.clusterPopoutModalData['selectedSmiles'] = undefined;
            this.clusterPopoutModalData['res'] = undefined;
            this.clusterPopoutModalData['group_id'] = undefined;
        },
        clusterPopoutModalIncGroupID: function() {
            var all_ids = this.clusteredResultsIndex[this.clusterPopoutModalData['selectedSmiles']];
            var idx = all_ids.indexOf(this.clusterPopoutModalData['group_id']);
            if (idx == all_ids.length-1) {
            } else {
                this.clusterPopoutModalData['group_id'] = all_ids[idx+1];
            }
            this.$forceUpdate();
        },
        clusterPopoutModalDecGroupID: function() {
            var all_ids = this.clusteredResultsIndex[this.clusterPopoutModalData['selectedSmiles']];
            var idx = all_ids.indexOf(this.clusterPopoutModalData['group_id']);
            if (idx == 0) {
            } else {
                this.clusterPopoutModalData['group_id'] = all_ids[idx-1];
            }
            this.$forceUpdate();
        },
        openClusterEditModal: function(selected, group_id) {
            if(selected == undefined) {
                alert('No target molecule selected. Please select a molecule in the tree.')
                return
            }
            /*
             * cannot deselect
            this.clearSelection();
            if (network) {
                network.unselectAll();
            }
            */
            if(group_id == undefined) {
                group_id = 0
            }
            this.$set(this.clusterEditModalData, 'selected', selected);
            this.$set(this.clusterEditModalData, 'selectedSmiles', selected.smiles);
            this.$set(this.clusterEditModalData, 'group_id', group_id);
            this.showClusterEditModal = true;
        },
        closeClusterEditModal: function() {
            this.showClusterEditModal = false;
            this.clusterEditModalData['selected'] = undefined;
            this.clusterEditModalData['selectedSmiles'] = undefined;
            this.clusterEditModalData['group_id'] = undefined;
        },
        clusteredit_dragstart_handler: function(precursor, event) {
            event.target.style.opacity = '0.4';
            event.dataTransfer.setData('text/plain', precursor.smiles);
            var img = new Image();
            img.src = this.getMolDrawEndPoint(precursor.smiles);
            // set opacity does not work..
            event.dataTransfer.setDragImage(img, 10, 10);
            event.dataTransfer.effectAllowed = 'all';
            // disable all buttons on dragging
            var buttons = document.querySelectorAll("button");
            buttons.forEach(function(e){e.style.pointerEvents = "none";});
        },
        clusteredit_dragend_handler: function(event) {
            event.target.style.opacity = '1';
            // enable all buttons
            var buttons = document.querySelectorAll("button");
            buttons.forEach(function(e){e.style.pointerEvents = "all";});
        },
        clusteredit_drop_handler: function(target, event) {
            event.preventDefault(); // important
            var s = event.dataTransfer.getData('text/plain'); // precursor.simles
            var r = this.results[this.clusterEditModalData['selectedSmiles']];
            // find precursor
            var old_gid;
            for (let x of r) {
                if (x.smiles == s) {
                    old_gid = x.group_id;
                    x.group_id = target.group_id;
                    break
                }
            }
            this.clusteredit_dragend_handler(event);
            clusteredit_dragleave_handler(event);
            this.detectClusterDeletion(this.clusterEditModalData['selectedSmiles'], old_gid);
        },
        clusteredit_drop_handler_newcluster: function(event) {
            event.preventDefault(); // important
            var s = event.dataTransfer.getData('text/plain'); // precursor.simles
            var r = this.results[this.clusterEditModalData['selectedSmiles']];
            var all_ids = this.clusteredResultsIndex[this.clusterEditModalData['selectedSmiles']];
            var new_gid = all_ids[all_ids.length-1]+1;
            var old_gid;
            for (let x of r) {
                if (x.smiles == s) {
                    old_gid = x.group_id;
                    x.group_id = new_gid;
                    break
                }
            }
            
            this.clusteredit_dragend_handler(event);
            clusteredit_dragleave_handler(event);
            this.detectClusterDeletion(this.clusterEditModalData['selectedSmiles'], old_gid);
        },
        detectClusterDeletion: function(selected, old_gid) {
            var all_ids = this.clusteredResultsIndex[selected];
            if (all_ids.indexOf(old_gid) == -1) {
                if (all_ids.length > 0) {
                    var idx = all_ids.findIndex(function(e){return e>old_gid});
                    if (idx == -1) idx = all_ids.length-1;
                    this.clusterEditModalData['group_id'] = all_ids[idx];
                } else {
                    this.clusterEditModalData['group_id'] = 0;
                }
                this.$forceUpdate();
            }
        },
        clusterEditModalIncGroupID: function() {
            var all_ids = this.clusteredResultsIndex[this.clusterEditModalData['selectedSmiles']];
            var idx = all_ids.indexOf(this.clusterEditModalData['group_id']);
            if (idx == all_ids.length-1) {
            } else {
                this.clusterEditModalData['group_id'] = all_ids[idx+1];
            }
            this.$forceUpdate();
        },
        clusterEditModalDecGroupID: function() {
            var all_ids = this.clusteredResultsIndex[this.clusterEditModalData['selectedSmiles']];
            var idx = all_ids.indexOf(this.clusterEditModalData['group_id']);
            if (idx == 0) {
            } else {
                this.clusterEditModalData['group_id'] = all_ids[idx-1];
            }
            this.$forceUpdate();
        },
        clusterEditModalDeletePrecursor: function(selected, smiles) {
            let res = confirm('This will remove the precursor completely and cannot be undone! Continue?')
            if (res) {
                var r = this.results[selected];
                var idx = r.findIndex(function(e){return e.smiles==smiles;});
                var old_gid = r[idx].group_id;
                r.splice(idx, 1);
                this.detectClusterDeletion(this.clusterEditModalData['selectedSmiles'], old_gid);
            }
        },
        // gid == undefined is to add a new cluster
        clusterEditModalAddPrecursor: function(selectedSmiles, smiles, gid) {
            var isshow = false;
            if (this.results[selectedSmiles] == undefined) {
                this.results[selectedSmiles] = [];
                gid = 0;
                isshow = true;
            }
            var all_ids = this.clusteredResultsIndex[selectedSmiles];
            if (gid == undefined) {
                isshow = true;
                if (all_ids.length == 0) {
                    gid = 0;
                } else {
                    gid = all_ids[all_ids.length-1]+1;
                }
            }
            var rank = 0;
            for (let i of this.results[selectedSmiles]) {
                rank = Math.max(rank, i.rank);
            }
            rank += 1;
            var r = {
                'show': isshow,
                'smiles': smiles,
                'smiles_split': smiles.split('.'),
                'group_id': gid,
                'score': undefined,
                'plausibility': undefined,
                'rank': rank,
                'num_examples': undefined,
                'necessary_reagent': undefined,
                'template_score': undefined,
                'templates': undefined,
            };
            this.results[selectedSmiles].push(r);
        },
        // if group_id == undefined, add to a new group
        openAddNewPrecursorModal: function(selectedSmiles, group_id) {
            this.showAddNewPrecursorModal = true;
            this.$set(this.addNewPrecursorModal, 'selectedSmiles', selectedSmiles == undefined ? this.selected.smiles : selectedSmiles);
            this.$set(this.addNewPrecursorModal, 'group_id', group_id == undefined ? 'undefined' : group_id.toString());
            this.$set(this.addNewPrecursorModal, 'newprecursorsmiles', '');
            this.$set(this.addNewPrecursorModal, 'nodupcheck', false);
        },
        closeAddNewPrecursorModal: function() {
            this.showAddNewPrecursorModal = false;
            this.addNewPrecursorModal['selectedSmiles'] = '';
            this.addNewPrecursorModal['group_id'] = '';
            this.addNewPrecursorModal['newprecursorsmiles'] = '';
            this.addNewPrecursorModal['nodupcheck'] = false;
        },
        checkDuplicatePrecursor: function(selectedSmiles, p) {
            var p_splited = new Set(p.split("."));
            for (s of this.results[selectedSmiles]) {
                var s_set = new Set(s['smiles_split']);
                if (subSet(s_set, p_splited) || subSet(p_splited, s_set)) {
                    return s;
                }
            }
            return undefined;
        },
        addNewPrecursorModalSubmit: async function() {
            var gid;
            if (this.addNewPrecursorModal['group_id'] == "undefined") {
                gid = undefined;
            } else {
                gid = Number(this.addNewPrecursorModal['group_id']);
            }
            try {
                isvalid = await this.validatesmiles(
                    this.addNewPrecursorModal['newprecursorsmiles'],
                    !this.allowResolve
                );
                if (!isvalid) {
                    this.addNewPrecursorModal['newprecursorsmiles'] =
                        await this.resolveChemName(
                            this.addNewPrecursorModal['newprecursorsmiles']
                        );
                }
            } catch(error) {
                var error_msg = 'unknown error';
                if ('message' in error) {
                    error_msg = error.name+':'+error.message;
                } else if (typeof(error) == 'string') {
                    error_msg = error;
                }
                alert('There was an error fetching precursors for this target with the supplied settings: '+error_msg);
                return
            }
            if (this.addNewPrecursorModal['newprecursorsmiles'] == undefined) {
                alert('There was an error during adding the precursor.');
            } else {
                if (!this.addNewPrecursorModal['nodupcheck']) {
                    var s = this.checkDuplicatePrecursor(
                        this.addNewPrecursorModal['selectedSmiles'],
                        this.addNewPrecursorModal['newprecursorsmiles']
                    );
                    if (s != undefined) {
                        alert('There may be a duplicated precursor: rank: '+s.rank+' cluster: '+s.group_id+'. If you still want to proceed, please select "No duplicate check" option.');
                        return
                    }
                }
                this.clusterEditModalAddPrecursor(
                    this.addNewPrecursorModal['selectedSmiles'],
                    this.addNewPrecursorModal['newprecursorsmiles'],
                    gid);
                this.$forceUpdate();
            }
        },
        getMolDrawEndPoint: function(precursor, isHighlight, isTransparent) {
            //  precursor can be
            //      1) a smiles string,
            //      2) a dict has properties "reacting_atoms" and "mapped_smiles"
            //      3) a dict has property "smiles"
            //  isTransparent is false by default
            //  isHighlight is set to this.isHighlight by default, but can be overidden
            if (isHighlight == undefined) {
                isHighlight = this.isHighlightAtom;
            }
            if (isTransparent == undefined) {
                isTransparent = false;
            }
            var smiles;
            var mapped_smiles;
            var reacting_atoms;
            if (typeof(precursor) == "string") {
                smiles = precursor;
                isHighlight = false;
            } else if (typeof(precursor) == "object") {
                if (precursor.mapped_smiles != undefined && precursor.reacting_atoms != undefined) {
                    mapped_smiles = precursor.mapped_smiles;
                    reacting_atoms = precursor.reacting_atoms;
                }
                if (precursor.smiles != undefined) {
                    smiles = precursor.smiles;
                }
            }
            if (isHighlight && mapped_smiles != undefined && reacting_atoms != undefined) {
                var res = '/draw/highlight/smiles='+encodeURIComponent(mapped_smiles)+'&reacting_atoms='+encodeURIComponent('['+reacting_atoms.toString()+']')+'&bonds=0'
            } else {
                if (smiles == undefined) {
                    console.log('Error: cannot plot precursor='+precursor)
                    return ''
                }
                var res = '/draw/smiles/' + encodeURIComponent(smiles);
            }
            if (isTransparent) {
                res += '?transparent=1';
            }
            return res;
        },
        isModalOpen: function() {
            var res = false;
            res = res || this.showSettingsModal;
            res = res || this.showDownloadModal;
            res = res || this.showLoadModal;
            res = res || this.showClusterPopoutModal;
            res = res || this.showClusterEditModal;
            res = res || this.showAddNewPrecursorModal;
            return res;
        },
        startTour: function() {
            if (this.target) {
                this.clear();
            }
            tour.init();
            tour.restart();
        },
        initClusterShowCard: function(selected) {
            // always sort first
            var reactionSorting = this.sortingCategory;
            // this.results[selected].sort(function(a, b) {
            //     var a_ = a[reactionSorting] == undefined ? 0 : a[reactionSorting];
            //     var b_ = b[reactionSorting] == undefined ? 0 : b[reactionSorting];
            //     return b_ - a_;
            // })
            // init show to false
            // init first reactionLimit clusters/precursors to true
            var numShow = 0;
            var visited_groups = new Set();
            for (precursor of this.results[selected]) {
                if (this.allowCluster) {
                    if (visited_groups.has(precursor.group_id)) {
                        this.$set(precursor, 'show', false);
                    } else {
                        this.$set(precursor, 'show', true);
                        visited_groups.add(precursor.group_id);
                    }
                } else { // !allowCluster
                    this.$set(precursor, 'show', true);
                }
            }
        },
        groupPrecursors: function(precursors) {
            var grouped = {};
            for (let i = 0; i < precursors.length; i++) {
                var precursor = precursors[i];
                if (grouped[precursor.group_id]) {
                    grouped[precursor.group_id].push(precursor);
                }
                else {
                    grouped[precursor.group_id] = new Array(precursor);
                }
            }
            return Object.values(grouped);
        },
        requestClusterId: function(selected) {
            showLoader();
            var all_smiles = [];
            var all_scores = [];
            var i;
            for (i = 0; i < this.results[selected].length; i++) {
                all_smiles.push(this.results[selected][i].smiles);
                var s = this.results[selected][i].score;
                if (s == undefined) {
                    s = 0;
                }
                all_scores.push(s);
            }
            var url = '/api/cluster/?';
            var params = {
                original:       selected,
                outcomes:       all_smiles,
                feature:        this.clusterOptions.feature,
                fp_name:        this.clusterOptions.fingerprint,
                fpradius:       this.clusterOptions.fpRadius,
                fpnbits:        this.clusterOptions.fpBits,
                cluster_method: this.clusterOptions.cluster_method,
                scores:          all_scores,
            };
            var queryString = Object.keys(params).map((key) => {
                return encodeURIComponent(key) + '=' + encodeURIComponent(params[key])
            }).join('&');
            
            fetch_param = {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/x-www-form-urlencoded',
                    'X-CSRFToken': getCookie('csrftoken'),
                },
                body: queryString,
            };
            
            fetch(url, fetch_param)
            .then(resp => {
                if (!resp.ok) {
                    throw resp.status;
                }
                return resp;
            })
            .then(resp => resp.json())
            .then(resp_json => {
                if ('error' in resp_json) {
                    hideLoader()
                    throw resp_json['error'];
                } else {
                    var group_ids = resp_json['group_id'];
                    var i;
                    for (i = 0; i < this.results[selected].length; i++) {
                        this.$set(this.results[selected][i], 'group_id', group_ids[i]);
                    }
                    hideLoader()
                }
            })
            .catch((error) => {
                hideLoader()
                var error_msg = 'unknown error'
                if (typeof(error) == 'number') {
                    error_msg = 'Error code: ' + error;
                } else if (typeof(error) == 'string') {
                    error_msg = error;
                } else if ('message' in error) {
                    error_msg = error.name+':'+error.message;
                }
                alert('There was an error fetching cluster results for this target with the supplied settings: '+error_msg)
            })
        },
        createTargetNode: function(target) {
            return {
                id: 0,
                smiles: target,
                image: this.getMolDrawEndPoint(target),
                shape: "image",
                borderWidth: 3,
                type: 'chemical',
                value: 15,
                mass: 2,
                color: {
                    border: '#000088'
                }
            }
        },
        alreadyAddedToresults: function(chemical, results) {
            for (var res of results) {
              if (chemical.reactant_smiles.join('.') == res.smiles) {
                  return true
              }
            }
            return false
        },
        addResultsFromTreeBuilder: function(graph, target) {
            this.target = target
            this.data.nodes = new vis.DataSet([])
            this.data.nodes.add(this.createTargetNode(target))
            this.data.edges = new vis.DataSet([]);
            for (smiles in graph) {
              if (!this.results[smiles]) {
                this.results[smiles] = new Array()
              }
              let rank = 1
              graph[smiles].sort((a, b) => b.template_score - a.template_score)
              for (chemical of graph[smiles]) {
                if (this.alreadyAddedToresults(chemical, this.results[smiles])) {
                    continue
                }
                this.results[smiles].push({
                  rank: rank,
                  smiles_split: chemical.reactant_smiles,
                  smiles: chemical.reactant_smiles.join('.'),
                  plausibility: chemical.plausibility,
                  visit_count: chemical.visit_count,
                  price: chemical.price,
                  estimate_price: chemical.estimate_price,
                  template_score: chemical.template_score,
                  score: chemical.template_score
                })
                rank += 1
              }
            }
        },
        addPathsFromTreeBuilder: function(trees) {
            for (tree of trees) {
                this.walkTree(tree, 0)
              }
            initializeNetwork(this.data, true);
            network.on('selectNode', this.showInfo);
            network.on('deselectNode', this.clearSelection);
            this.toggleHierarchical()
        },
        walkTree: function(obj, parent) {
            var node = undefined
            for (graphNode of this.data.nodes.get()) {
                if (graphNode.graphId == obj.id) {
                    let graphParent = parentOf(graphNode.id, this.data.nodes, this.data.edges)
                    if (graphParent == -1) {
                        node = graphNode
                    } else {
                        graphParent = this.data.nodes.get(graphParent)
                        if ((graphParent.graphId == parent.graphId) && (graphParent.id == parent.id)) {
                            node = graphNode
                            break
                        }
                    }
                }
                if (graphNode.id == 0 && graphNode.smiles == obj.smiles) {
                    node = graphNode
                    continue
                }
            }
            if (!node) {
                node = ctaToNode(obj, this.data.nodes.length)
                if (node.type == 'reaction') {
                    for (templateId of node.templateIds) {
                        this.apiTemplateCount(templateId)
                    }
                    smi_split = node.reactionSmiles.split('>>')
                    for (reaction of this.results[smi_split[1]]) {
                        if (reaction.smiles == smi_split[0]) {
                            node.rank = reaction.rank
                            node.label = '#'+node.rank
                            reaction.inViz = true
                            break
                        }
                    }
                }
                this.data.nodes.add(node)
                if (parent != 0) {
                    this.data.edges.add({
                        from: parent.id,
                        to: node.id,
                        length: 0
                    })
                }
            }
            for (child of obj.children) {
                this.walkTree(child, node)
            }
        },
        loadFromTreeBuilder: function(objectId, numTrees) {
            this.allowCluster = false
            showLoader()
            fetch('/api/get-result/?id='+objectId)
            .then(resp => resp.json())
            .then(json => {
              if (json.error) {
                  alert(json.error)
              }
              var result = json['result'];
              var target = result['settings']['smiles'];
              var stats = result['result']['status'];
              var trees = result['result']['paths'];
              var graph = result['result']['graph'];
              this.addResultsFromTreeBuilder(graph, target)
              if (numTrees == 'all') {
                this.addPathsFromTreeBuilder(trees)
              }
              else {
                this.addPathsFromTreeBuilder(trees.slice(0, Number(numTrees)))
              }
            })
            .finally(() => hideLoader())
        }
    },
    computed: {
        // {'target_smiles0':[[{result0}, {result1}, ...], [...]], ...}
        clusteredResults: function() {
            var res = {};
            var x;
            for (x in this.results) {
                res[x] = this.groupPrecursors(this.results[x]);
                vueApp = this
                res[x].sort(function(a, b) {
                    let maxPropA = Math.max.apply(Math, a.map(function(obj) {
                        return obj[vueApp.sortingCategory]
                    }))
                    let maxPropB = Math.max.apply(Math, b.map(function(obj) {
                        return obj[vueApp.sortingCategory]
                    }))
                    return maxPropB - maxPropA
                })
            }
            return res;
        },
        // {'target_smiles0':[all possible unique group_ids sorted in accending order], ...}
        clusteredResultsIndex: function() {
            var res = {};
            var x;
            for (x in this.results) {
                var ids = new Set();
                for (let i of this.results[x]) {
                    ids.add(i.group_id);
                }
                res[x] = Array.from(ids).sort(function(a, b){return a-b});
            }
            return res;
        }
    },
    delimiters: ['%%', '%%'],
});

var tour = new Tour({
    storage: false,
    steps: [
        {
            title: "A guided tour through retrosynthesis",
            content: "Welcome to this guided tour through retrosynthesis planning using our interactive path planning tool. This will demonstrate the purpose of the tool and explain the user interface using a real example. Thanks to <a href='http://bootstraptour.com/' target='_blank'>bootstrap-tour</a> for the great guided tour JavaScript package making it very easy to provide this tour to you!",
            orphan: true,
            backdropContainer: '#body'
        },
        {
            element: "#target",
            title: "Start with a target compound",
            content: "You can start the retrosynthetic planning with a target compound and typing it's SMILES formatted string here. If the name resolver is enabled (see server icon to the right; click icon to toggle), you can also enter a chemical name. The name will be resolved using a third-party server (PubChem). For this tutorial we're going to explore an example reaction for <a href='https://en.wikipedia.org/wiki/Fluconazole' target='_blank'>Fluconazole</a>. Press 'Next' to continue!",
            placement: "bottom",
            onNext: function() {
                app.target = 'OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F'
            }
        },
        {
            element: "#target",
            title: "Fluconazole",
            content: "Here's the SMILES string for Fluconazole. If you're unfamiliar with the SMILES format, try using a software like ChemDraw to draw a structure and copy it's SMILES string (right click -> molecule -> copy as -> SMILES). Click next to continue!",
            placement: "bottom",
            onNext: function() {
                if (app.data.nodes.length == null | app.data.nodes.length == 0) {
                    app.changeTarget();
                }
            }
        },
        {
            element: "#network",
            title: "One-step retrosynthesis results",
            content: "When the results are ready, they will be shown in the main window. The target molecule you entereted will be shown in the middle inside a <span class='blue-text'>blue</span> box (it is currently selected). You can click and drag on empty space in the window to translate the entire network. You can try to rearrange nodes by clicking and dragging on them. Take a second to enjoy the inverted gravity model, courtesy of <a href='http://visjs.org' target='_blank'>vis.js</a>. Scrolling inside the window will zoom in and out.",
            placement: 'right',
            backdropContainer: '#network'
        },
        {
            element: "#network",
            title: "Predicted reactions",
            content: "The children nodes of your target molecule (one is highlighted, for example) represent predicted <b>reactions</b> that may result in your target molecule. The number inside this node is the rank of the precursor, scored by the precursor prioritization method currently selected (more on this later).",
            onShown: function () {
                network.selectNodes([1]);
                app.selected = app.data.nodes.get(1);
            },
            placement: 'right',
        },
        {
            element: "#network",
            title: "Reactants",
            content: "The children node(s) of <b>reactions</b> represent <b>chemicals</b>, and are the predicted reactants for this reaction. Chemicals in a <span class='red-text'>red</span> box weren't found in our buyables database. <b>Chemicals</b> in a <span class='green-text'>green</span> box are buyable.",
            placement: 'right',
            onNext: function() {
                app.data.nodes.forEach(function(n) {
                    if (n.smiles == 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1') {
                        network.selectNodes([n.id])
                        app.selected = app.data.nodes.get(n.id);
                    }
                })
            }
        },
        {
            element: "#network",
            placement: 'right',
            title: "Reactants",
            content: "For this example, we'll see if we can predict a reaction to make this reaction's non-buyable reactant (it's been selected for you) from buyable starting materials."
        },
        {
            element: '#expand-btn',
            title: "Expanding chemical nodes",
            content: "The non-buyable chemical for which we'd like to make a new prediction has been highlighted for you. Next you'd click the <b>Expand Node</b> button. Click next to see what happens when you click this button.",
            placement: "bottom",
            reflex: true,
            onNext: function() {
                app.expandNode();
            }
        },
        {
            element: '#network',
            title: "Expanding chemical nodes",
            content: "A new prediction was made for this non-buyable chemical, and when everything is ready the results will be added to the network visualization. It might look hectic at first, but the appropriate node positions should resolve quickly (thanks again to the <a href='http://visjs.org' target='_blank'>vis.js</a> inverted gravity!). If not, click and drag a node to give it a jiggle.",
            placement: "right"
        },
        {
            element: '#details',
            title: "Result details",
            content: "You may have noticed there's been a lot going on on the right side of the screen in addition to the changes in the graph visualization. On this side, details of the currently selected node are shown. In this case, a <b>chemical</b> node is selected. At the top you can see its SMILES string, its cost in $/g and a 2d rendering of its structure.",
            placement: "left"
        },
        {
            element: '#details',
            title: "Precursors",
            content: "Additionally, if you've already made a retrosynthetic prediction for the currently selected <b>chemical</b>, you'll see list of the precursor results. Each entry shows the reactants for the reaction to make the currently selected chemical with some additional information such as a relative score and the number of examples there were for the templates that support the suggested reaction. You can reorder these results by each metric using the drop-down menu above. If you haven't performed a retrosynthetic prediction for the selected chemical, the same <b>Expand Node</b> button you used before will be shown.",
            placement: "left"
        },
        {
            element: '#details',
            title: "Adding and removing reactions",
            content: "You may also notice there are many more precursor results shown on the right side here than were added into the graph visualization (it's a scrolling list) - this is to keep things tidy in the visualization. By default, only the top 5 results (scored by retro'score') are added to the visualization (this can be changed in the settings menu). The plus (+) and minus (-) buttons can be used to add and remove each reaction to the visualization. Go ahead and give it a try if you'd like.",
            placement: "left",
            onNext: function() {
                app.data.nodes.forEach(function(n) {
                    if (n.reactionSmiles == 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1.c1nc[nH]n1>>OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F') {
                        network.selectNodes([n.id])
                        app.selected = app.data.nodes.get(n.id);
                    }
                })
            }
        },
        {
            element: '#details',
            title: "Viewing reaction details",
            content: "If you have a reaction node selected, the right side of your screen will show you details for that reaction. At the top you can see the reaction SMILES, a 2d rendering, and similar reaction scores that you have seen before. You will also see a list of links to templates that support the reaction. Clicking one will open a new tab with more details about each template. There is also a link to 'Evaluate reaction in new tab', which will let you predict reaction conditions and evaluate the reaction in the forward direction.",
            placement: "left",
            onNext: function() {
                app.data.nodes.forEach(function(n) {
                    if (n.smiles == 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1') {
                        network.selectNodes([n.id])
                        app.selected = app.data.nodes.get(n.id);
                    }
                })
            }
        },
        {
            element: "#network",
            title: "Understanding the network",
            content: "We can see that the prediction gave us a few reactions that use buyable starting materials (<b>reaction</b> nodes with children <b>chemical</b> nodes highlighted in green) to make the new target we were interested in. Now we have a full path to our original target, Fluconazole, starting with buyable compounds.",
            placement: "right"
        },
        {
            element: "#expand-btn",
            title: "Other buttons",
            content: "In addition to expanding nodes, you can easily delete selected nodes, or children of a selected node using the corresponding red buttons above the graph visualization. The 'Toggle cluster' button will group the currently selected node and its children into one node cluster (this may be useful to keep things organized). Clicking this button with a cluster selected will expand it to show all of the nodes again.",
            placement: "bottom"
        },
        {
            element: "#settings-btn",
            title: "Getting more/less results",
            content: "By default, only the top 5 predicted reactions for each target are shown. If you want more or less, open this settings window."
        },
        {
            element: "#download-btn",
            title: "Saving results",
            content: "You can save the network structure (JSON of nodes and edges) and download it to your computer."
        },
        {
            element: "#load-btn",
            title: "Restoring results",
            content: "You can restore a previously saved network here."
        },
        {
            title: "End of tour",
            content: "That's the end of the guided tour. Go ahead and change the target, and build your own reaction networks, or continue to expand this one further.",
            orphan: true
        }
    ]
});

function closeAll() {
    app.showSettingsModal = false;
    app.showLoadModal = false;
    app.showDownloadModal = false;
    app.showClusterPopoutModal = false;
    app.showClusterEditModal = false;
    app.showAddNewPrecursorModal = false;
}

/* key binding */
var keys = vis.keycharm();
keys.bind("esc", closeAll, 'keyup');
