function showLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "block";
}

function hideLoader() {
    var loader = document.getElementsByClassName("loader")[0];
    loader.style.display = "none";
}

var app = new Vue({
    el: '#app',
    data: {queues: []},
    created: function() {
        this.update();
    },
    methods: {
        update: function() {
            showLoader();
            fetch('/api/celery/')
            .then(resp => resp.json())
            .then(json => {
                console.log(json['queues']);
                this.queues = json['queues'];
            })
            .finally(() => hideLoader())
        }
    },
    delimiters: ['%%', '%%'],
});