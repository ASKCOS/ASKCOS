var container = document.getElementsByClassName('container')[0];
container.style.width=null;

function showLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "block";
}

function hideLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "none";
}

function addReactions(reactions, sourceNode, nodes, edges, reactionLimit) {
    var reactionSorting = "retroscore";
    reactions.sort(function(a, b) {
        return b[reactionSorting] - a[reactionSorting]
    })
    for (n in reactions) {
        console.log(n)
        if (n >= reactionLimit) {
            break;
        }
        addReaction(reactions[n], sourceNode, nodes, edges)
    }
}

function addReaction(reaction, sourceNode, nodes, edges) {
    var rId = nodes.max('id').id+1;
    nodes.add({
        id: rId,
        label: '#'+reaction['rank'],
        rank: reaction['rank'],
        ffScore: reaction['plausibility'].toFixed(3),
        retroscore: reaction['score'].toFixed(3),
        templateScore: reaction['template_score'].toFixed(3),
        numExamples: reaction['num_examples'],
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
        fetch('/ajax/price_smiles/?smiles='+encodeURIComponent(smi))
        .then(resp => resp.json())
        .then(json => {
            var mysmi = json['smiles'];
            var ppg = json['ppg'];
            if (ppg == 0) {
                ppg = "not buyable"
            }
            var buyable = (json['ppg']!=0);
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
                image: window.location.origin+"/draw/smiles/"+encodeURIComponent(mysmi),
                shape: "image",
                borderWidth: 2,
                type: 'chemical',
                mass: 1,
                value: 10,
                ppg: ppg,
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

function initializeNetwork(data) {
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
        }
    };
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

Vue.component('modal', {
    template: '#modal-template'
})

var app = new Vue({
    el: '#app',
    data: {
        target: '',
        network: {},
        data: {
            nodes: {},
            edges: {}
        },
        results: {},
        templateNumExamples: {},
        nodeStructure: {},
        allowResolve: false,
        showSettingsModal: false,
        showLoadModal: false,
        showDownloadModal: false,
        downloadName: "network.json",
        modalData: {},
        selected: null,
        reactionLimit: 5,
        templatePrioritization: "Relevance",
        precursorScoring: "RelevanceHeuristic",
        numTemplates: 1000,
        maxCumProb: 0.999,
        minPlausibility: 0.01,
        sortingCategory: "score"
    },
    beforeMount: function() {
        this.allowResolve = this.$el.querySelector('[ref="allowResolve"]').checked;
    },
    methods: {
        requestUrl: function(smiles) {
            var url = '/api/retro/?';
            var params = {
                target: smiles,
                template_prioritization: this.templatePrioritization,
                precursor_prioritization: this.precursorScoring,
                num_templates: this.numTemplates,
                max_cum_prob: this.maxCumProb,
                filter_threshold: this.minPlausibility,
            }
            var queryString = Object.keys(params).map((key) => {
                return encodeURIComponent(key) + '=' + encodeURIComponent(params[key])
            }).join('&');
            return url+queryString;
        },
        resolveTarget: function() {
            if (this.allowResolve) {
                var url = 'https://cactus.nci.nih.gov/chemical/structure/'+encodeURIComponent(this.target)+'/smiles'
                console.log(url)
                fetch(url)
                    .then(resp => {
                        if (resp.status == 404) {
                            throw Error(resp.statusText);
                        }
                        else {
                            return resp.text()
                        }
                    })
                    .then(text => {
                        console.log(text);
                        this.target = text;
                        this.changeTarget();
                    })
                    .catch(err => {
                        alert('Cannot resolve "'+this.target+'" to smiles');
                })
            }
            else {
                this.changeTarget();
            }
        },
        changeTarget: function() {
            showLoader();
            var url = this.requestUrl(this.target);
            fetch(url)
                .then(resp => {
                    if (!resp.ok) {
                        throw Error(resp.statusText);
                    }
                    return resp;
                })
                .then(resp => resp.json())
                .then(json => {
                    this.data.nodes = new vis.DataSet([
                        {
                            id: 0,
                            smiles: this.target,
                            image: window.location.origin+"/draw/smiles/"+encodeURIComponent(this.target),
                            shape: "image",
                            borderWidth: 3,
                            type: 'chemical',
                            value: 15,
                            mass: 2,
                            color: {
                                border: '#000088'
                            }
                        }
                    ])
                    this.data.edges = new vis.DataSet([]);
                    initializeNetwork(this.data);
                    network.on('selectNode', this.showInfo);
                    network.on('deselectNode', this.clearSelection);
                    this.results[this.target] = json['precursors'];
                    addReactions(json['precursors'], this.data.nodes.get(0), this.data.nodes, this.data.edges, this.reactionLimit);
                    this.getTemplateNumExamples(json['precursors']);
                    hideLoader();
                    fetch('/api/price/?smiles='+encodeURIComponent(this.target))
                        .then(resp => resp.json())
                        .then(json => {
                            var ppg = json['price'];
                            if (ppg == 0) {
                                ppg = "not buyable"
                            }
                            this.data.nodes.update({id: 0, ppg: ppg});
                            network.selectNodes([0]);
                            this.selected = this.data.nodes.get(0);
                    })
                })
                .catch(error => {
                    hideLoader();
                    alert('There was an error fetching precursors for this target with the supplied settings')
                })
        },
        expandNode: function() {
            if (this.isModalOpen() || typeof(network) == "undefined") {
                return
            }
            showLoader();
            var selected = network.getSelectedNodes();
            console.log(selected.length);
            if (selected.length == 0) {
              hideLoader();
            }
            for (n in selected) {
                var nodeId = selected[n];
                var node = this.data.nodes.get(nodeId)
                if (node.type != 'chemical') {
                    alert('Cannot expand reaction; try expanding with a chemical node selected');
                    continue
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
                        this.results[smi] = reactions;
                        if (reactions.length==0) {
                            alert('No precursors found!')
                        }
                        addReactions(reactions, this.data.nodes.get(nodeId), this.data.nodes, this.data.edges, this.reactionLimit);
                        this.getTemplateNumExamples(reactions);
                        this.selected = node;
                        this.reorderResults();
                        hideLoader();
                    })
                    .catch(error => {
                        hideLoader();
                        alert('There was an error fetching precursors for this target with the supplied settings')
                    })
            }
        },
        getTemplateNumExamples: function(reactions) {
            console.log(reactions)
            for (reaction of reactions) {
                for (templateId of reaction['templates']) {
                    if (typeof(this.templateNumExamples[templateId]) == 'undefined') {
                        fetch('/api/template/?id='+templateId)
                        .then(resp => resp.json())
                        .then(json => {
                            var id = json["request"]["id"][0];
                            var count = json["template"]["count"];
                            this.templateNumExamples[id] = count;
                        })
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
                    alert('We do not reccomend deleting chemical nodes! It will leaving its parent reaction node missing information! Only delete reaction nodes with this button.')
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
                    alert('We do not reccomend deleting children of reaction nodes! It will leaving the reaction node missing information! Only delete children of chemical nodes with this button.')
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
        },
        clearSelection: function() {
            this.selected = null;
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
            results.sort((a, b) => b[sortingCategory] - a[sortingCategory])
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
            this.clearSelection();
            if (network) {
                network.unselectAll();
            }
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
        isModalOpen: function() {
            return app.showSettingsModal || app.showDownloadModal || app.showLoadModal
        },
        startTour: function() {
            this.clear();
            tour.init();
            tour.restart();
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
            content: "You can start the retrosynthetic planning with a target compound and typing it's SMILES formatted string here. If the name resolver is enabled (see server icon to the right; click icon to toggle), you can also enter a chemical name. The name will be resolved using a third-party server (NIH CACTUS). For this tutorial we're going to explore an example reaction for <a href='https://en.wikipedia.org/wiki/Fluconazole' target='_blank'>Fluconazole</a>. Press 'Next' to continue!",
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
}

var keys = vis.keycharm();
keys.bind("esc", closeAll, 'keyup');
