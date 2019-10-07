var container = document.getElementsByClassName('container')[0];
container.style.width=null;

var loader = document.getElementsByClassName("loader")[0];
loader.style.top = "300px";

function showLoader() {
    loader.style.display = "block";
}
showLoader();

function hideLoader() {
    loader.style.display = "none";
}

function hideNetwork(n) {
  var networkDiv = document.querySelectorAll('.tree-graph')[n];
  var hideDiv = document.querySelectorAll('.hideNetwork')[n];
  var showDiv = document.querySelectorAll('.showNetwork')[n];
  networkDiv.style.display = 'none';
  hideDiv.style.display = 'none';
  showDiv.style.display = '';
}

function hideAllNetworks() {
  for (var n=0; n<document.querySelectorAll('.tree-graph').length; n++) {
    hideNetwork(n)
  }
}

function showNetwork(n) {
  var networkDiv = document.querySelectorAll('.tree-graph')[n];
  var hideDiv = document.querySelectorAll('.hideNetwork')[n];
  var showDiv = document.querySelectorAll('.showNetwork')[n];
  networkDiv.style.display = '';
  hideDiv.style.display = '';
  showDiv.style.display = 'none';
}

function showAllNetworks() {
  for (var n=0; n<document.querySelectorAll('.tree-graph').length; n++) {
    showNetwork(n)
  }
}

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

function colorOf(child) {
  if (child['ppg']) {
    if (child['as_reactant'] || child['as_product']) {
      return "#1B5E20" // green
    }
    else {
      return '#FFC400' // yellow
    }
  }
  else {
    if (child['as_reactant'] || child['as_product']) {
      return  '#E65100' // orange
    }
    else {
      return '#B71C1C' // red
    }
  }
}

function makeNode(child, id) {
  var node = {};
  if (child['is_chemical']) {
    node['id'] = id
    node['ppg'] = child['ppg']
    node['smiles'] = child['smiles']
    node['image'] = "/draw/smiles/"+child['smiles']
    node['shape'] = "image"
    node['type'] = 'chemical'
    var buyableString = (Number(child['ppg'])) ? '$'+child['ppg']+'/g' : 'not buyable'
    node['as_reactant'] = child['as_reactant']
    node['as_product'] = child['as_product']
    node['title'] = `${child['smiles']}<br>
    ${child['as_reactant']} precedents as reactant<br>
    ${child['as_product']} precedents as product<br>
    ${buyableString}`
    node['borderWidth'] = 2
    node['color'] = {
      border: colorOf(child)
    }
  }
  else if (child['is_reaction']) {
    node['type'] = 'reaction'
    node['id'] = id
    node['necessary_reagent'] = child['necessary_reagent']
    node['num_examples'] = child['num_examples']
    node['plausibility'] = child['plausibility']
    node['smiles'] = child['smiles']
    node['template_score'] = child['template_score']
    node['tforms'] = child['tforms']
    node['label'] = `${child['num_examples']} examples
FF score: ${Number(child['plausibility']).toFixed(3)}`
    node['font'] = {align: 'center'}
  }
  return node
}

function makeTreeData(obj, parent, nodes, edges) {
  var node = makeNode(obj, nodes.length);
  nodes.add(node);
  if (parent != 0) {
    edges.add({
      from: parent.id,
      to: node.id,
      length: 0
    })
  }
  for (child of obj.children) {
    makeTreeData(child, node, nodes, edges);
  }
}

function makeNodes(children, root) {
  var nodes = [];
  var edges = [];
  var childrenOfChildren = [];
  for (child of children) {
    nodes.push(makeNode(child))
    edges.push({
      from: root.id,
      to: child.id
    })
    childrenOfChildren.push(child.children)
  }
  var result = {
    nodes: nodes,
    edges: edges,
    children: childrenOfChildren
  }
  return result
}

function initializeNetwork(data, elementDiv) {
    var container = elementDiv;
    var options = {
        nodes: {
            color: {
                border: '#000000',
                background: '#FFFFFF'
            },
            shapeProperties: {
                useBorderWithImage: true,
                useImageSize: true
            }
        },
        edges: {
            length: 1
        },
        interaction: {
            multiselect: false,
            hover: true,
            dragNodes: false,
            dragView: false,
            selectConnectedEdges: false,
            tooltipDelay: 0,
            zoomView: false
        },
        layout: {
          hierarchical: {
            levelSeparation: 150,
            nodeSpacing: 175,
          }
        },
        physics: false
    };
    network = new vis.Network(container, data, options);
    network.on("beforeDrawing",  function(ctx) {
        ctx.save();
        ctx.setTransform(1, 0, 0, 1, 0, 0);
        ctx.fillStyle = '#ffffff';
        ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height)
        ctx.restore();
    })
    network.storePositions();
    var nodes = data.nodes;
    var first = nodes.min('y');
    var last = nodes.max('y');
    var y = last.y - first.y;
    container.style.height = `${y+200}px`;
    var first = nodes.min('x');
    var last = nodes.max('x');
    var x = last.x - first.x;
    container.style.width = `${x+200}px`;
    return network
}

function getCookie(name) {
    var cookieValue = null;
    if (document.cookie && document.cookie !== '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = cookies[i].trim();
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}
var csrftoken = getCookie('csrftoken');

var app = new Vue({
    el: '#app',
    data: {
      networks: [],
      networkData: [],
      resultId: "",
      numChemicals: 0,
      numReactions: 0,
      trees: [],
      settings: {},
      selected: null
    },
    mounted: function() {
      this.resultId = this.$el.getAttribute('data-id');
      this.getResult(this.resultId);
    },
    methods: {
      getResult: function(id) {
        showLoader();
        fetch('/api/get-result/?id='+id)
          .then(resp => resp.json())
          .then(json => {
            var result = json['result'];
            var stats = result['result']['status'];
            var trees = result['result']['paths'];
            this.numChemicals = stats[0];
            this.numReactions = stats[1];
            this.trees = trees;
            this.settings = result['settings'];
            this.renderSettings();
            this.buildTrees(this.trees);
            sleep(1000).then(() => app.fitNetworks()).then(() => hideLoader());
          })
      },
      buildTrees: function(trees) {
        var container = document.getElementById('left-pane')
        var statsDiv = document.createElement('div')
        statsDiv.classList.add('text-center')
        statsDiv.style.margin = "10px auto"
        statsDiv.innerHTML = `After expanding (with ${this.settings.known_bad_reactions.length} banned reactions, ${this.settings.forbidden_molecules.length} banned chemicals), ${this.numChemicals} total chemicals and ${this.numReactions} total reactions`;
        container.appendChild(statsDiv);
        var hideShowAllDiv = document.createElement('div');
        hideShowAllDiv.classList.add('text-center');
        hideShowAllDiv.style.margin = "10px auto";
        hideShowAllDiv.innerHTML = '<a onclick=hideAllNetworks()>Hide all</a>';
        hideShowAllDiv.innerHTML += ' ';
        hideShowAllDiv.innerHTML += '<a onclick=showAllNetworks()>Show all</a>';
        container.appendChild(hideShowAllDiv);
        var count = 0
        for (tree of trees) {
          count += 1;
          var titleDiv = document.createElement('div');
          titleDiv.classList.add('text-center');
          titleDiv.innerHTML = '<h3 style="display: inline; margin-right: 5px">Option '+count+'</h3>';
          titleDiv.innerHTML += '<a class="hideNetwork" onclick=hideNetwork('+(count-1)+')>hide</a>'
          titleDiv.innerHTML += '<a style="display: none" class="showNetwork" onclick=showNetwork('+(count-1)+')>show</a>'
          container.appendChild(titleDiv);
          var divElem = document.createElement('div');
          divElem.classList.add('tree-graph');
          container.appendChild(divElem);
          this.buildTree(tree, divElem, count);
        }
      },
      buildTree: function(tree, elem, id) {
        var nodes = new vis.DataSet([]);
        var edges = new vis.DataSet([]);
        makeTreeData(tree, 0, nodes, edges);
        var data = {
          nodes: nodes,
          edges: edges
        }
        var network = initializeNetwork(data, elem);
        network.id = id;
        network.on('selectNode', function(params) {
          app.showNode(id, params.nodes[0])
        });
        network.fit();
        this.networks.push(network);
        this.networkData.push(data);
      },
      resizeCanvases: function() {
        for (var n in this.networks) {
          var network = this.networks[n];
          network.storePositions();
          var nodes = this.networkData[n].nodes;
          var first = nodes.min('y');
          var last = nodes.max('y');
          var x = last.x - first.x;
          var y = last.y - first.y;
          var div = document.querySelectorAll('.tree-graph')[n];
          div.style.height = `${y+200}px`;
        }
      },
      fitNetworks: function() {
        for (var network of this.networks) {
          network.fit();
        }
      },
      renderSettings: function() {
        console.log(this.settings);
        var settingsDiv = document.createElement('div');
        var targetDiv = document.createElement('div');
        targetDiv.innerHTML = '<div class="text-center">Target: '+this.settings.smiles+'</div><div class="text-center"><img src="/draw/smiles/'+this.settings.smiles+'"></div>'
        settingsDiv.appendChild(targetDiv);
        var settingsTitle = document.createElement('h3');
        settingsTitle.innerHTML = 'Settings:';
        settingsDiv.appendChild(settingsTitle);
        var settingsTable = document.createElement('table');
        settingsTable.classList.add('table');
        settingsTable.innerHTML = '<tr><th>Expansion settings:</th><td>Max. depth: '+this.settings.max_depth+'</td><td>Max. branching factor: '+this.settings.max_branching+'</td></tr>'
        settingsDiv.appendChild(settingsTable);
        settingsTable.innerHTML += '<tr><th></th><td>Num. templates: '+this.settings.template_count+'</td><td>Max cum. prob: '+this.settings.max_cum_template_prob+'</td></tr>'
        settingsTable.innerHTML += '<tr><th></th><td colspan=2>Expansion time (s): '+this.settings.expansion_time+'</td></tr>'
        settingsTable.innerHTML += '<tr><th>Stop criteria:</th><td colspan=2>Maximum chemical price ($/g): '+this.settings.max_ppg+'</td></tr>'
        if (this.settings.max_natom_dict.logic != 'none') {
          settingsTable.innerHTML += `<tr><th></th><td colspan=2>Chemical property logic: C=${this.settings.max_natom_dict.C} N=${this.settings.max_natom_dict.N} H=${this.settings.max_natom_dict.H} O=${this.settings.max_natom_dict.O}</td></tr>`
        }

        if (this.settings.min_chemical_history_dict.logic != 'none') {
          settingsTable.innerHTML += `<tr><th></th><td colspan=2>Chemical popularity logic: Min. freq. as reactant=${this.settings.min_chemical_history_dict.as_reactant} Min. freq. as product=${this.settings.min_chemical_history_dict.as_product}</td></tr>`
        }
        settingsTable.innerHTML += `<tr>
          <th>Evaluation settings:</th>
          <td>Min. plausibility: ${this.settings.filter_threshold}</td>
          <td></td>
        </tr>`
        document.getElementById('settings').appendChild(settingsDiv);
      },
      blacklist: function() {
        n = 0
        for (network of this.networks) {
          var nodeId = network.getSelectedNodes();
          if (nodeId.length) {
            var desc = prompt("Please enter a reason (for your records only)", "no reason");
            var now = Date.now();
            var datetime = now.toString('MMMM dd, yyyy, hh:mm:ss tt');
            var nodes = this.networkData[n].nodes;
            node = nodes.get(nodeId[0]);
            if (node['type'] == 'chemical') {
              var url =   '/ajax/user_blacklist_chemical/'
            }
            else {
              var url = '/ajax/user_blacklist_reaction/'
            }
            $.ajax({
                type: 'POST',
                url: url,
                data: {
                    smiles: node.smiles,
                    csrfmiddlewaretoken: csrftoken,
                    desc: desc,
                    datetime: datetime,
                },
                dataType: 'json',
                success: function (data) {
                    if (data.err) {
                        alert(data.err);
                    } else {
                        alert('Blacklisted chemical "' + node.smiles + '" at ' + datetime);
                    }
                }
            })
            break
          }
          n += 1;
        }
      },
      showNode: function(networkId, nodeId) {
        for (network of this.networks) {
          if (network.id == networkId) {
            nodeId = network.getSelectedNodes()[0]
            this.selected = this.networkData[networkId-1].nodes.get(nodeId)
          }
          else {
            network.unselectAll()
          }
        }
      }
    },
    delimiters: ['%%', '%%'],
});
