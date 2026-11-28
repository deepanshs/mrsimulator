/* The basicstrap theme defines DOCUMENTATION_OPTIONS inline, in the form of Sphinx 1.
 * The search and navigation scripts of Sphinx 6 and later read further options, and
 * the root of the documentation from the html element. Both are set here. */
(function () {
  var options = window.DOCUMENTATION_OPTIONS || (window.DOCUMENTATION_OPTIONS = {});
  var defaults = {
    LANGUAGE: document.documentElement.lang || "en",
    BUILDER: "html",
    LINK_SUFFIX: options.FILE_SUFFIX || ".html",
    SOURCELINK_SUFFIX: ".txt",
    NAVIGATION_WITH_KEYS: false,
    SHOW_SEARCH_SUMMARY: true,
    ENABLE_SEARCH_SHORTCUTS: true,
  };
  Object.keys(defaults).forEach(function (key) {
    if (options[key] === undefined) options[key] = defaults[key];
  });
  if (document.documentElement.dataset.content_root === undefined) {
    document.documentElement.dataset.content_root = options.URL_ROOT || "./";
  }
})();
