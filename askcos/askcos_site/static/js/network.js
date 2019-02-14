function showLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "block";
}

function hideLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "none";
}

function addReactions(reactions, sourceNode, nodes, edges, reactionLimit, reactionSorting) {
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
        numExamples: reaction['num_examples'],
        templateIds: reaction['tforms'],
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
    for (smi of reaction['smiles_split']) {
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
                image: window.location.origin+"/draw/smiles/"+mysmi,
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
            children.push(e.to)
            for (child of allChildrenOf(e.to, nodes, edges)) {
                children.push(child)
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
        nodeStructure: {},
        showReactionModal: false,
        showChemicalModal: false,
        showSettingsModal: false,
        showLoadModal: false,
        showDownloadModal: false,
        downloadName: "network.json",
        modalData: {},
        reactionLimit: 5,
        reactionSorting: "retroscore",
        templatePrioritization: "Relevance",
        precursorScoring: "RelevanceHeuristic",
        numTemplates: 100,
        minPlausibility: 0.75
    },
    methods: {
        requestUrl: function(smiles) {
            var url = '/api/retro/?';
            var params = {
                smiles: smiles,
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
                            image: window.location.origin+"/draw/smiles/"+this.target,
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
                    addReactions(json['precursors'], this.data.nodes.get(0), this.data.nodes, this.data.edges, this.reactionLimit, this.reactionSorting);
                    hideLoader();
                })
                .catch(error => {
                    hideLoader();
                    alert('There was an error fetching precursors for this target with the supplied settings')
                })
        },
        expandNode: function() {
            showLoader();
            var selected = network.getSelectedNodes();
            for (nodeId of selected) {
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
                        if (reactions.length==0) {
                            alert('No precursors found!')
                        }
                        addReactions(reactions, this.data.nodes.get(nodeId), this.data.nodes, this.data.edges, this.reactionLimit, this.reactionSorting);
                        hideLoader();
                    })
                    .catch(error => {
                        hideLoader();
                        alert('There was an error fetching precursors for this target with the supplied settings')
                    })
            }
        },
        deleteNode: function() {
            var selected = network.getSelectedNodes();
            for (nodeId of selected) {
                removeChildrenFrom(nodeId, this.data.nodes, this.data.edges);
                this.data.nodes.remove(nodeId);
            }
            cleanUpEdges(this.data.nodes, this.data.edges);
        },
        deleteChildren: function() {
            var selected = network.getSelectedNodes();
            for (nodeId of selected) {
                removeChildrenFrom(nodeId, this.data.nodes, this.data.edges);
            }
            cleanUpEdges(this.data.nodes, this.data.edges);
        },
        nodeInfo: function() {
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
            for (childId of childrenId) {
                reactants.push(this.data.nodes.get(childId).smiles)
            }
            var rxnSmiles = reactants.join('.')+'>>'+product
            this.modalData['product'] = product;
            this.modalData['reactionSmiles'] = rxnSmiles;
            this.modalData['reactionSmilesImg'] = "/draw/reaction/"+rxnSmiles;
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
            }})(file);
            reader.readAsText(file)
        },
        clear: function() {
            this.target = '';
            this.data.nodes.remove(this.data.nodes.getIds())
            this.data.edges.remove(this.data.edges.getIds())
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
            content: "Welcome to this guided tour through retrosynthesis planning using our reaction network explorer. This will deomnstrate the purpose of the tool and explain the user interface using a real example. Thanks to <a href='http://bootstraptour.com/' target='_blank'>bootstrap-tour</a> for the great guided tour JavaScript package making it very easy to provide this tour to you!",
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
            content: "Here's the SMILES string for Fluconazole. If you're unfamiliar with the SMILES format, try using a software like ChemDraw to draw a structure and get it's SMILES string. Click the 'Submit' button then click next to continue!",
            placement: "bottom",
            onNext: function() {
                if (app.data.nodes.length == null | app.data.nodes.length == 0) {
                    app.changeTarget();
                }
            }
        },
        {
            title: "One-step retrosynthesis results",
            content: "When the results are ready, they will be shown in the main window. The target molecule you entereted will be shown in the middle inside a <span class='blue-text'>blue</span> box. You can click and drag on empty space in the window to translate the entire network. You can try to rearrange nodes by clicking and dragging on them. Take a second to enjoy the inverted gravity model, courtesy of <a href='http://visjs.org' target='_blank'>vis.js</a>. Scrolling inside the window will zoom in and out.",
            orphan: true,
            backdropContainer: '#network'
        },
        {
            title: "Predicted reactions",
            content: "The children nodes of your target molecule (one is highlighted, for example) represent predicted <b>reactions</b> that result in your target molecule. The number inside this node represents a fast filter score related to reaction likelyhood.",
            onShown: function () {
                network.selectNodes([1])
            },
            orphan: true
        },
        {
            title: "Reactants",
            content: "The children node of the <b>reaction</b> represent <b>chemicals</b>, and are the predicted reactants for this reaction. Chemicals in a <span class='red-text'>red</span> box weren't found in our buyables database. Chemicals in a <span class='green-text'>green</span> box are buyable. For this example, we'll see if we can predict a reaction to make this reaction's non-buyable chemical from buyable starting materials.",
            orphan: true,
            onNext: function() {
                app.data.nodes.forEach(function(n) {
                    if (n.smiles == 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1') {
                        network.selectNodes([n.id])
                    }
                })
            }
        },
        {
            element: '#expand-btn',
            title: "Expanding chemical nodes",
            content: "The non-buyable chemical for which we'd like to make a new prediction has been highlighted for you. Click this 'Expand node' button to run a one-step retrosynthetic prediction for the selected chemical and add the results to your reaction network.",
            placement: "bottom",
            reflex: true
        },
        {
            element: '#network',
            title: "Expanding chemical nodes",
            content: "When the results are ready they'll be added to the network. It might look hectic at first, but appropriate node positions should resolve quickly. If not, click and drag a node to give it a jiggle.",
            placement: "top"
        },
        {
            title: "Understanding the network",
            content: "We can see that the prediction gave us a few reactions that use buyable starting materials (<b>reaction</b> nodes with children <b>chemical</b> nodes highlighted in green) to make the new target we were interested in. Now we have a full path to our original target, Fluconazole, starting with buyable compounds.",
            orphan: true
        },
        {
            title: "Other buttons",
            content: "Some of the other buttons are self-explanatory. You can delete selected nodes, delete children of selected nodes, collapse/cluster nodes and their children, and view more detailed reaction information.",
            orphan: true
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

var keys = vis.keycharm();
keys.bind("backspace", app.deleteNode, 'keyup');
keys.bind("delete", app.deleteNode, 'keyup');