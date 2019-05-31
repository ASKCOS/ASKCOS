function showPopup() {
  document.querySelector('.support-popup').style.display = 'block';
}

function hidePopup() {
  document.querySelector('.support-popup').style.display = 'none';
}

function submitSupport() {
  var mailtoString = "mailto:incoming+mlpds-mit-askcos-askcos-12564933-issue-@incoming.gitlab.com"
  mailtoString += "?subject="+document.querySelector('#support-module').value
  mailtoString += "|"+document.querySelector('#support-category').value
  if (document.querySelector('#support-shared').checked) {
    mailtoString += "|shared"
  }
  else {
    mailtoString += "|"
  }
  mailtoString += "| "+document.querySelector('#support-subject').value
  mailtoString += "&body=Please add your comments here. Please don't modify the beginning of the subject line; we use this for internal purposes"
  window.location = mailtoString;
}

document.querySelector('#support-anchor').onclick = showPopup;
document.querySelector('#support-submit').onclick = submitSupport;
document.querySelector('#support-cancel').onclick = hidePopup;
