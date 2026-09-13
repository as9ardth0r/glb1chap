async function main() {
  const grid = document.getElementById("grid");
  const empty = document.getElementById("empty");
  const meta = document.getElementById("meta");
  const sortSelect = document.getElementById("sort-by");
  const novelOnly = document.getElementById("novel-only");

  let payload;
  try {
    const res = await fetch("data/molecules.json");
    if (!res.ok) throw new Error("HTTP " + res.status);
    payload = await res.json();
  } catch (err) {
    meta.textContent = "Impossible de charger data/molecules.json";
    empty.style.display = "block";
    return;
  }

  if (!payload.molecules || payload.molecules.length === 0) {
    meta.textContent = "Cible : " + payload.target + " — 0 candidat";
    empty.style.display = "block";
    return;
  }

  meta.textContent = `Cible : ${payload.target} — ${payload.count} candidat(s) conforme(s) au TPP — généré le ${new Date(payload.generated_at).toLocaleString("fr-FR")}`;

  function render() {
    const key = sortSelect.value;
    let molecules = [...payload.molecules];

    if (novelOnly.checked) {
      molecules = molecules.filter(m => m.is_novel === true);
    }

    molecules.sort((a, b) => {
      const av = a[key], bv = b[key];
      if (av === null || av === undefined) return 1;
      if (bv === null || bv === undefined) return -1;
      // docking_score : plus négatif = meilleur, donc tri croissant
      if (key === "docking_score") return av - bv;
      return bv - av;
    });

    grid.innerHTML = "";
    if (molecules.length === 0) {
      empty.style.display = "block";
      return;
    }
    empty.style.display = "none";

    for (const m of molecules) {
      const card = document.createElement("div");
      card.className = "card";

      const rGroup = m.recipe ? m.recipe.r_group : "?";
      const dockingText = m.docking_score !== null && m.docking_score !== undefined
        ? `${m.docking_score} kcal/mol` : "non docké";
      const noveltyText = m.is_novel === true ? "nouveau"
        : m.is_novel === false ? `connu${m.chemical_name ? " (" + m.chemical_name + ")" : ""}`
        : "non vérifié";

      card.innerHTML = `
        <div class="svg-wrap">${m.svg || "<em>pas de dépiction</em>"}</div>
        <span class="r-group">N-${rGroup}</span>
        <h3>${m.id}</h3>
        <div class="stats">
          <div>MW: <b>${m.mw}</b></div>
          <div>LogP: <b>${m.logp}</b></div>
          <div>TPSA: <b>${m.tpsa}</b></div>
          <div>QED: <b>${m.qed}</b></div>
          <div>Docking: <b>${dockingText}</b></div>
          <div>Nouveauté: <b>${noveltyText}</b></div>
        </div>
        <div class="fitness-bar" style="width:${Math.round((m.fitness || 0) * 100)}%"></div>
      `;
      grid.appendChild(card);
    }
  }

  sortSelect.addEventListener("change", render);
  novelOnly.addEventListener("change", render);
  render();
}

main();
