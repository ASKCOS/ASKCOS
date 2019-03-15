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
        label: reaction['plausibility'].toFixed(3),
        ffScore: reaction['plausibility'].toFixed(3),
        retroscore: reaction['score'].toFixed(3),
        templateScore: reaction['template_score'].toFixed(3),
        numExamples: reaction['num_examples'],
        templateIds: reaction['templates'],
        reactionSmiles: reaction.smiles+'>>'+sourceNode.smiles,
        type: 'reaction',
        value: 1,
        mass: 0.5
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
        color: '#000000'
    })
    for (n in reaction['smiles_split']) {
        var smi = reaction['smiles_split'][n];
        fetch('/ajax/price_smiles/?smiles='+encodeURIComponent(smi))
        .then(resp => resp.json())
        .then(json => {
            var mysmi = json['smiles'];
            var ppg = json['ppg'];
            if (ppg == 0) {
                ppg = "N/A"
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
                type: 'chemical',
                mass: 0.5,
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
                color: '#000000'
            })
        })
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
        // save current translate/zoom
        ctx.save();
        // reset transform to identity
        ctx.setTransform(1, 0, 0, 1, 0, 0);
        // fill background with solid white
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height)
        // restore old transform
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
        nodeStructure: {},
        allowResolve: false,
        showSettingsModal: false,
        showLoadModal: false,
        showDownloadModal: false,
        showKeyboardOverlay: false,
        downloadName: "network.json",
        modalData: {},
        selected: null,
        reactionLimit: 5,
        templatePrioritization: "Relevance",
        precursorScoring: "RelevanceHeuristic",
        numTemplates: 100,
        minPlausibility: 0.75
    },
    beforeMount: function() {
        this.allowResolve = this.$el.querySelector('[ref="allowResolve"]').value;
    }
    methods: {
        requestUrl: function(smiles) {
            var url = '/api/retro/?';
            var params = {
                target: smiles,
                template_prioritization: this.templatePrioritization,
                precursor_prioritization: this.precursorScoring,
                num_templates: this.numTemplates,
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
                fetch(url)
                    .then(resp => resp.text())
                    .then(text => {
                        this.target = text;
                        this.changeTarget();
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
                            type: 'chemical',
                            value: 15,
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
                    hideLoader();
                    fetch('/api/price/?smiles='+encodeURIComponent(this.target))
                        .then(resp => resp.json())
                        .then(json => {
                            var ppg = json['price'];
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
            if (this.isModalOpen()) {
                return
            }
            showLoader();
            var selected = network.getSelectedNodes();
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
                        this.selected = node;
                        hideLoader();
                    })
                    .catch(error => {
                        hideLoader();
                        alert('There was an error fetching precursors for this target with the supplied settings')
                    })
            }
        },
        deleteNode: function() {
            if (this.isModalOpen()) {
                return
            }
            var selected = network.getSelectedNodes();
            for (n in selected) {
                var nodeId = selected[n];
                removeChildrenFrom(nodeId, this.data.nodes, this.data.edges);
                this.data.nodes.remove(nodeId);
            }
            cleanUpEdges(this.data.nodes, this.data.edges);
        },
        deleteChildren: function() {
            var selected = network.getSelectedNodes();
            for (n in selected) {
                var nodeId = selected[n];
                removeChildrenFrom(nodeId, this.data.nodes, this.data.edges);
            }
            cleanUpEdges(this.data.nodes, this.data.edges);
        },
        nodeInfo: function() {
            if (this.isModalOpen()) {
                return
            }
            var selected = network.getSelectedNodes();
            var node = this.data.nodes.get(selected[0]);
            if (node.type == 'reaction') {
                this.reactionInfo(node);
            }
            else if (node.type == 'chemical') {
                this.chemicalInfo(node);
            }
            else {
                alert('Cannot show information for this node.')
            }
        },
        reactionInfo: function(node) {
            var reactants = [];
            var childrenId = childrenOf(node.id, this.data.nodes, this.data.edges);
            var parentId = parentOf(node.id, this.data.nodes, this.data.edges);
            var product = this.data.nodes.get(parentId).smiles
            for (n in childrenId) {
                var childId = childrenId[n];
                reactants.push(this.data.nodes.get(childId).smiles)
            }
            var rxnSmiles = reactants.join('.')+'>>'+product
            this.modalData['product'] = product;
            this.modalData['reactionSmiles'] = rxnSmiles;
            this.modalData['reactionSmilesImg'] = "/draw/reaction/"+encodeURIComponent(rxnSmiles);
            this.modalData['node'] = node;
            this.showReactionModal = true;
        },
        chemicalInfo: function(node) {
            this.modalData['node'] = node;
            this.showChemicalModal = true;
        },
        download: function() {
            if (this.data.nodes.length == null) {
                alert("There's no network to download!")
                return
            }
            var downloadData = {nodes: [], edges: []}
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
            addReaction(reaction, selected, this.data.nodes, this.data.edges)
        },
        remFromResults: function(selected, reaction) {
            var rsmi = reaction.smiles+'>>'+selected.smiles;
            var selectedChildren = this.data.nodes.get(childrenOf(selected.id, this.data.nodes, this.data.edges));
            for (var child of selectedChildren) {
                console.log(child, rsmi);
                if (child.reactionSmiles == rsmi) {
                    removeChildrenFrom(child.id, this.data.nodes, this.data.edges);
                    this.data.nodes.remove(child.id);
                    cleanUpEdges(this.data.nodes, this.data.edges);
                    break;
                }
            }
        },
        showInfo: function(obj) {
            var nodeId = obj.nodes[obj.nodes.length-1];
            var node = this.data.nodes.get(nodeId);
            this.selected = node;
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
            else if (modalName == "keyboard") {
                this.showKeyboardOverlay = true
            }
        },
        isModalOpen: function() {
            return app.showSettingsModal || app.showChemicalModal || app.showDownloadModal || app.showLoadModal || app.showReactionModal || app.showKeyboardOverlay
        },
        startTour: function() {
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
            content: "Welcome to this guided tour through retrosynthesis planning using our reaction network explorer. This will demonstrate the purpose of the tool and explain the user interface using a real example. Thanks to <a href='http://bootstraptour.com/' target='_blank'>bootstrap-tour</a> for the great guided tour JavaScript package making it very easy to provide this tour to you!",
            orphan: true,
            backdropContainer: '#body'
        },
        {
            element: "#target",
            title: "Start with a target compound",
            content: "You can start the retrosynthetic planning with a target compound and typing it's SMILES formatted string here. For this tutorial we're going to explore an example reaction for <a href='https://en.wikipedia.org/wiki/Fluconazole' target='_blank'>Fluconazole</a>. Press 'Next' to continue!",
            placement: "top",
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
            content: "The children nodes of your target molecule (one is highlighted, for example) represent predicted <b>reactions</b> that may result in your target molecule. The number inside this node represents a fast filter score, or plausibility, related to reaction likelihood.",
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
            content: "You may have noticed there's been a lot going on on the right side of the screen in addition to the changes in the graph visualization. On this side, details of the currently selected node are shown. In this case, a <b>chemical</b> node is selected. At the top you can see its SMILES string, its cost in $/g and a 2d rendering of its structure. Note: a price of N/A means it was not found in our buyables database.",
            placement: "left"
        },
        {
            element: '#details',
            title: "Precursors",
            content: "Additionally, if you've already made a retrosynthetic prediction for the currently selected <b>chemical</b>, you'll see list of the precursor results. Each entry shows the reactants for the reaction to make the currently selected chemical with some additional information such as a relative score and the number of examples there were for the templates that support the suggested reaction. If you haven't performed a retrosynthetic prediction for the selected chemical, the same <b>Expand Node</b> button you used before will be shown.",
            placement: "left"
        },
        {
            element: '#details',
            title: "Adding and removing reactions",
            content: "You may also notice there are many more precursor results shown on the right side here than were added into the graph visualization (it's a scrolling list) - this is to keep things tidy in the visualization. By default, only the top 5 results are added to the visualization (this can be changed in the settings menu). The plus (+) and minus (-) buttons can be used to add and remove each reaction to the visualization. Go ahead and give it a try if you'd like.",
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
            placement: "top"
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
            element: "#keyboard-overlay",
            title: "Keyboard shortcuts",
            content: "There are a few keyboard shortcuts to help expand or delete nodes that are explained here."
        },
        {
            title: "End of tour",
            content: "That's the end of the guided tour. Go ahead and change the target, and build your own reaction networks, or continue to expand this one further.",
            orphan: true
        }
    ]
});

function closeAll() {
    app.showReactionModal = false;
    app.showChemicalModal = false;
    app.showSettingsModal = false;
    app.showLoadModal = false;
    app.showDownloadModal = false;
    app.showKeyboardOverlay = false;
}

var keys = vis.keycharm();
keys.bind("esc", closeAll, 'keyup');
keys.bind("backspace", app.deleteNode, 'keyup');
keys.bind("delete", app.deleteNode, 'keyup');
keys.bind("d", app.deleteNode, 'keyup');
keys.bind("e", app.expandNode, 'keyup');
keys.bind("i", app.nodeInfo, 'keyup');