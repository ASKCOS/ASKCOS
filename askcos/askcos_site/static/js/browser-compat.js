function unsupportedAlert() {
    alert('Your browser does not support features this application requires! Specifically, the application is written using ES6 JavaScript features. Please use a recent version of Chrome, Firefox, Edge, Safari, or Opera.')
}

function checkFeatures() {
    if (typeof(fetch) == 'undefined') {
        unsupportedAlert();
        return
    }
    try {new Function("(x) => (x)")}
    catch(e) {
        unsupportedAlert();
        return
    }
}
checkFeatures();