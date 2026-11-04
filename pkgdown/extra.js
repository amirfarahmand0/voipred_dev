// The "On this page" menu of pkgdown (bootstrap-toc) shows two heading levels.
// Add the third (e.g. 4.1.1 under 4.1): nested like the second level, shown
// when its parent section is active, and followed by the scrollspy.
$(function () {
  var $toc = $("#toc");
  if (!$toc.length) return;
  var parent = null;
  $("main").find("h2, h3, h4").not("[data-toc-skip]").each(function () {
    var tag = this.tagName.toLowerCase();
    if (tag !== "h4") {
      parent = tag === "h3" ? this : null;
      return;
    }
    if (!parent || !parent.id) return;
    var $parentLink = $toc.find('a.nav-link[href="#' + parent.id + '"]');
    if (!$parentLink.length) return;
    if (!this.id) this.id = "section-" + Math.random().toString(36).slice(2, 10);
    var $list = $parentLink.next("ul.nav");
    if (!$list.length) {
      $list = $('<ul class="nav navbar-nav"></ul>');
      $parentLink.after($list);
    }
    var text = $(this).clone().find(".anchor").remove().end().text().trim();
    $list.append($("<li></li>").append(
      $('<a class="nav-link"></a>').attr("href", "#" + this.id).text(text)));
  });
  if (window.bootstrap && bootstrap.ScrollSpy) {
    var spy = bootstrap.ScrollSpy.getInstance(document.body);
    if (spy) spy.refresh();
  }
});
