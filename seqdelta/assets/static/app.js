document.addEventListener("DOMContentLoaded", () => {
  const filterRoot = document.querySelector("[data-filter-root]");
  if (filterRoot) {
    const typeFilter = filterRoot.querySelector("[data-filter-type]");
    const effectFilter = filterRoot.querySelector("[data-filter-effect]");
    const searchFilter = filterRoot.querySelector("[data-filter-search]");
    const table = document.querySelector("[data-filter-table]");
    const rows = Array.from(table.querySelectorAll("tbody tr"));

    const applyFilters = () => {
      const typeValue = typeFilter.value.trim().toLowerCase();
      const effectValue = effectFilter.value.trim().toLowerCase();
      const searchValue = searchFilter.value.trim().toLowerCase();

      rows.forEach((row) => {
        const rowType = (row.dataset.type || "").toLowerCase();
        const rowEffect = (row.dataset.effect || "").toLowerCase();
        const rowText = row.textContent.toLowerCase();
        const visible =
          (!typeValue || rowType === typeValue) &&
          (!effectValue || rowEffect === effectValue) &&
          (!searchValue || rowText.includes(searchValue));
        row.style.display = visible ? "" : "none";
      });
    };

    typeFilter.addEventListener("change", applyFilters);
    effectFilter.addEventListener("change", applyFilters);
    searchFilter.addEventListener("input", applyFilters);
  }

  const triggers = Array.from(document.querySelectorAll("[data-view-trigger]"));
  const panels = Array.from(document.querySelectorAll("[data-view-panel]"));
  if (triggers.length && panels.length) {
    triggers.forEach((trigger) => {
      trigger.addEventListener("click", () => {
        const target = trigger.dataset.viewTrigger;
        triggers.forEach((button) => button.classList.toggle("is-active", button === trigger));
        panels.forEach((panel) => panel.classList.toggle("is-hidden", panel.dataset.viewPanel !== target));
      });
    });
  }
});
