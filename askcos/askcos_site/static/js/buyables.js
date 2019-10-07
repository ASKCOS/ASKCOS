function showLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "block";
}

function hideLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "none";
}

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


Vue.component('modal', {
    template: '#modal-template'
})

var app = new Vue({
    el: '#app',
    data: {
        buyables: [],
        uploadFile: '',
        searchSmilesQuery: '',
        searchSourceQuery: '',
        searchExact: true,
        searchLimit: 100,
        showAddModal: false,
        showUploadModal: false,
        addBuyableSmiles: '',
        addBuyablePrice: 1,
        addBuyableSource: '',
        uploadFileFormat: 'json'
    },
    methods: {
        search: function() {
            showLoader()
            fetch('/api/buyables/search?q='+encodeURIComponent(this.searchSmilesQuery)+'&source='+this.searchSourceQuery+'&exact='+this.searchExact+'&limit='+this.searchLimit)
            .then(resp => resp.json())
            .then(json => {
                console.log(json['buyables']);
                this.buyables = json['buyables'];
            })
            .finally(() => hideLoader())
        },
        handleFileUpload: function() {
            this.uploadFile = this.$refs.file.files[0]
            if (this.uploadFile.name.endsWith('.json')) {
                this.uploadFileFormat = 'json'
            }
            if (this.uploadFile.name.endsWith('.csv')) {
                this.uploadFileFormat = 'csv'
            }
        },
        handleUploadSubmit: function() {
            showLoader()
            if (this.uploadFile == '') {
                alert('Please select a file to upload')
                hideLoader()
                return
            }
            let formData = new FormData()
            formData.append('file', this.uploadFile)
            formData.append('format', this.uploadFileFormat)
            fetch('/api/buyables/upload', 
                {
                    'method': 'POST',
                    headers: {
                        'X-CSRFToken': getCookie('csrftoken')
                    },
                    body: formData
                }
            )
            .then(resp => resp.json())
            .then(json => {
                if (json.error) {
                    alert(json.error)
                    hideLoader()
                    return
                }
                alert('Successfully added/modified '+json.count+' of '+json.total+' entries')
                if (json.added.length > 0) {
                    this.buyables.unshift(...json.added)
                }
                if (json.updated.length > 0) {
                    for (updated of json.updated) {
                        let inList = false
                        for (buyable of this.buyables) {
                            if (buyable._id == updated._id) {
                                inList = true
                                buyable.ppg = updated.ppg
                                buyable.source = updated.source
                                break
                            }
                        }
                        if (!inList) {
                            this.buyables.unshift(updated)
                        }
                    }
                }
            })
            .finally(() => {
                this.uploadFile = ''
                hideLoader()
            })
        },
        addBuyable: function() {
            showLoader()
            fetch('/api/buyables/add?smiles='+encodeURIComponent(this.addBuyableSmiles)+'&ppg='+this.addBuyablePrice+'&source='+this.addBuyableSource)
            .then(resp => resp.json())
            .then(json => {
                if (json.error || json.success != true) {
                    alert('Error adding buyable compound')
                }
                else {
                    if (json.inserted) {
                        this.buyables.unshift(json.inserted)
                    }
                    if (json.updated) {
                        for (buyable of this.buyables) {
                            if (buyable._id == json.updated._id) {
                                buyable.ppg = json.updated.ppg
                                buyable.source = json.updated.source
                            }
                        }
                    }
                }
            })
            .finally(() => hideLoader())
        },
        deleteBuyable: function(_id) {
            showLoader()
            fetch('/api/buyables/delete?_id='+encodeURIComponent(_id))
            .then(resp => resp.json())
            .then(json => {
                if (json.error) {
                    alert(json.error)
                }
                if (json['success']) {
                    for (i in this.buyables) {
                        if (this.buyables[i]['_id']==_id) {
                            this.buyables.splice(i, 1)
                        }
                    }
                }
            })
            .finally(() => hideLoader())
        }
    },
    delimiters: ['%%', '%%'],
});