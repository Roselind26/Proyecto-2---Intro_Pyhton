# -*- coding: utf-8 -*-
"""
GC% en 16S rRNA — Dashboard interactivo (Proyecto 2)

Funcionalidades principales:
- Cargar archivos FASTA de genes 16S rRNA.
- Calcular %GC y %N de cada secuencia.
- Aplicar criterios de calidad (longitud y porcentaje de N).
- Resumir resultados por organismo.
- Generar un gráfico de barras HORIZONTAL comparando el %GC.
- Generar conclusiones automáticas:
    • Conclusión general (reglas simples).
    • Conclusión generada por un agente de IA usando OpenAI (gpt-4o-mini).

Requisitos de instalación (una sola vez):
    py -3.12 -m pip install streamlit pandas matplotlib openai

Para ejecutar el dashboard:
    py -3.12 -m streamlit run streamlit_app.py

La API key de OpenAI se debe guardar en la variable de entorno:
    OPENAI_API_KEY
"""

import io
import os
import re
from typing import List, Dict, Tuple

import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
from openai import OpenAI  # Cliente oficial de OpenAI

# -------------------------------------------------------------------
# CONFIGURACIÓN GENERAL DEL DASHBOARD
# -------------------------------------------------------------------
st.set_page_config(page_title="GC% 16S rRNA - Dashboard", layout="wide")
plt.rcParams["font.family"] = "Times New Roman"

st.title("GC% en 16S rRNA — Dashboard interactivo")
st.caption(
    "Sube archivos FASTA, ajusta los umbrales de calidad, compara el %GC por organismo "
    "y genera conclusiones automáticas (reglas + agente IA con OpenAI)."
)

# -------------------------------------------------------------------
# CONFIGURACIÓN DEL CLIENTE OPENAI (AGENTE IA)
# -------------------------------------------------------------------
def get_openai_client() -> OpenAI | None:
    """
    Crea un cliente de OpenAI usando la variable de entorno OPENAI_API_KEY.
    Si no está configurada o algo falla, devuelve None.
    """
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        return None
    try:
        client = OpenAI(api_key=api_key)
        return client
    except Exception:
        return None


def generar_conclusion_openai(df_sum: pd.DataFrame) -> str:
    """
    Usa el modelo gpt-4o-mini de OpenAI para generar una conclusión en español
    basada en la tabla de resumen por organismo.

    Si no hay API key o la petición falla, devuelve un mensaje explicativo.
    """
    client = get_openai_client()
    if client is None:
        return (
            "La funcionalidad de IA (OpenAI) no está disponible porque no se encontró "
            "una API key válida en la variable de entorno 'OPENAI_API_KEY'."
        )

    if df_sum.empty or df_sum["GC_promedio"].isna().all():
        return "La IA no puede generar una conclusión porque no hay datos válidos en el resumen por organismo."

    # Seleccionamos columnas relevantes para dar contexto a la IA
    cols = [c for c in df_sum.columns if c in ["Organismo", "GC_promedio", "N_promedio", "Longitud_total_ACGT"]]
    tabla_markdown = df_sum[cols].to_markdown(index=False)

    prompt = f"""
Eres un bioinformático experto en análisis de 16S rRNA y contenido GC.

A continuación tienes una tabla de resumen por organismo, donde:
- GC_promedio = porcentaje promedio de GC en el gen 16S rRNA.
- N_promedio  = porcentaje promedio de bases ambiguas N.
- Longitud_total_ACGT = suma de longitudes efectivas de las secuencias analizadas.

Tabla de datos:

{tabla_markdown}

Con base en estos resultados, escribe una conclusión en ESPAÑOL que:
1) Indique claramente qué organismo tiene el mayor %GC y cuál el menor.
2) Describa el rango de variación del %GC entre los organismos.
3) Comente brevemente qué implicaciones biológicas generales tiene un %GC más alto o más bajo
   (estabilidad del ADN/ARN, genomas reducidos, nichos ecológicos, etc.).
4) Use un tono adecuado para el cierre de un informe académico corto.
5) Tenga entre 6 y 10 líneas de texto en un solo párrafo (no usar viñetas).

No inventes organismos que no aparezcan en la tabla.
"""

    try:
        response = client.chat.completions.create(
            model="gpt-4o-mini",  # modelo ligero y rápido
            messages=[
                {
                    "role": "system",
                    "content": "Eres un asistente experto en bioinformática y redacción científica."
                },
                {"role": "user", "content": prompt},
            ],
            temperature=0.4,
            max_tokens=600,
        )
        # Nuevo SDK: el contenido está en message.content
        texto = response.choices[0].message.content
        return texto
    except Exception as e:
        return f"No se pudo generar la conclusión con OpenAI (IA): {e}"

# -------------------------------------------------------------------
# CONJUNTO DE ORGANISMOS DEL PROYECTO
# -------------------------------------------------------------------

# 9 organismos base del enunciado
BASE_9 = {
    "Mycobacterium tuberculosis",
    "Corynebacterium diphtheriae",
    "Bacillus subtilis",
    "Escherichia coli",
    "Salmonella enterica",
    "Pseudomonas aeruginosa",
    "Clostridium botulinum",
    "Borrelia burgdorferi",
    "Mycoplasma genitalium",
}

# 5 organismos adicionales sugeridos (para mayor exigencia)
EXTRA_5 = {
    "Staphylococcus aureus",
    "Vibrio cholerae",
    "Helicobacter pylori",
    "Streptococcus pneumoniae",
    "Listeria monocytogenes",
}

# Tabla de accesiones NCBI para 16S rRNA
ACCESIONES: Dict[str, str] = {
    "Mycobacterium tuberculosis": "NR_102810.2",
    "Corynebacterium diphtheriae": "NR_037079.1",
    "Bacillus subtilis": "NR_112116.2",
    "Escherichia coli": "NR_024570.1",
    "Salmonella enterica": "NR_074910.1",
    "Pseudomonas aeruginosa": "NR_074828.1",
    "Clostridium botulinum": "NR_029157.1",
    "Borrelia burgdorferi": "NR_044732.2",
    "Mycoplasma genitalium": "NR_026155.1",
    # 5 adicionales:
    "Staphylococcus aureus": "NR_113956.1",
    "Vibrio cholerae": "NR_118568.1",
    "Helicobacter pylori": "NR_028911.1",
    "Streptococcus pneumoniae": "NR_114926.1",
    "Listeria monocytogenes": "NR_074996.1",
}

# -------------------------------------------------------------------
# FUNCIONES DE UTILIDAD (LECTURA FASTA Y CÁLCULOS)
# -------------------------------------------------------------------
def read_fasta_text(text: str) -> str:
    """
    Lee contenido FASTA en texto plano:
    - Elimina líneas de encabezado (empiezan con '>').
    - Elimina caracteres que NO sean A/C/G/T/N (mayúscula o minúscula).
    - Devuelve una única secuencia concatenada en mayúsculas.
    """
    seq_lines: List[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            continue  # saltar encabezados
        clean = re.sub(r"[^ACGTNacgtn]", "", line.strip())
        seq_lines.append(clean)
    return "".join(seq_lines).upper()


def gc_n_content(seq: str) -> Tuple[float, float, int, int, int, int, int, int]:
    """
    Calcula:
    - %GC sobre bases A/C/G/T.
    - %N sobre todas las bases (A+C+G+T+N).
    Devuelve:
      (%GC, %N, A, C, G, T, N_total, length_acgt)
    """
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    n = seq.count("N")

    acgt = a + c + g + t
    total = acgt + n

    gc_pct = float("nan") if acgt == 0 else (g + c) * 100.0 / acgt
    n_pct = float("nan") if total == 0 else n * 100.0 / total

    return gc_pct, n_pct, a, c, g, t, n, acgt


def infer_organism_from_filename(name: str) -> str:
    """
    Intenta inferir el nombre del organismo a partir del nombre de archivo.
    Ejemplo:
        'Escherichia_coli_16S.fasta' → 'Escherichia coli'
    """
    base = name.rsplit(".", 1)[0]
    base = re.sub(r"[_\-]+", " ", base).strip()
    toks = base.split()
    if len(toks) >= 2:
        return toks[0].capitalize() + " " + toks[1].lower()
    return base


def build_df_from_uploads(files, len_min: int, len_max: int, amb_threshold: float) -> pd.DataFrame:
    """
    Construye un DataFrame a nivel de ENTRADA (archivo FASTA).
    Para cada archivo:
      - Lee la secuencia.
      - Calcula %GC y %N.
      - Aplica flags de longitud y ambigüedad.
    """
    rows = []
    for uf in files:
        text = uf.getvalue().decode("utf-8", errors="ignore")
        seq = read_fasta_text(text)
        gc_pct, n_pct, a, c, g, t, n, length = gc_n_content(seq)
        organismo = infer_organism_from_filename(uf.name)
        acceso = ACCESIONES.get(organismo, "")

        # Flags de calidad
        flag_len = not (len_min <= length <= len_max)
        flag_amb = (n_pct >= amb_threshold) if pd.notna(n_pct) else True

        rows.append({
            "Organismo": organismo,
            "Acceso_NCBI": acceso,
            "Archivo": uf.name,
            "Longitud_bases_ACGT": length,
            "A": a,
            "C": c,
            "G": g,
            "T": t,
            "N": n,
            "%GC": round(gc_pct, 3) if pd.notna(gc_pct) else gc_pct,
            "%N": round(n_pct, 3) if pd.notna(n_pct) else n_pct,
            "Flag_Longitud": flag_len,
            "Flag_Ambiguedad": flag_amb,
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("%GC", ascending=False).reset_index(drop=True)
    return df


def summarize_by_organism(df_entries: pd.DataFrame) -> pd.DataFrame:
    """
    Construye el resumen a nivel de ORGANISMO:
      - n_entradas
      - Longitud_total_ACGT
      - GC_promedio
      - N_promedio
      - Flags agregados (any)
    """
    if df_entries.empty:
        return pd.DataFrame(columns=[
            "Organismo", "n_entradas", "Longitud_total_ACGT",
            "GC_promedio", "N_promedio",
            "Flag_Longitud_any", "Flag_Ambiguedad_any",
        ])

    df_sum = (
        df_entries.groupby("Organismo", as_index=False)
        .agg({
            "Archivo": "count",
            "Longitud_bases_ACGT": "sum",
            "%GC": "mean",
            "%N": "mean",
            "Flag_Longitud": "any",
            "Flag_Ambiguedad": "any",
        })
        .rename(columns={
            "Archivo": "n_entradas",
            "Longitud_bases_ACGT": "Longitud_total_ACGT",
            "%GC": "GC_promedio",
            "%N": "N_promedio",
            "Flag_Longitud": "Flag_Longitud_any",
            "Flag_Ambiguedad": "Flag_Ambiguedad_any",
        })
        .sort_values("GC_promedio", ascending=False, ignore_index=True)
    )

    df_sum["GC_promedio"] = df_sum["GC_promedio"].round(3)
    df_sum["N_promedio"] = df_sum["N_promedio"].round(3)

    return df_sum

# -------------------------------------------------------------------
# GRÁFICO: BARRAS HORIZONTALES (MEJOR PARA NOMBRES LARGOS)
# -------------------------------------------------------------------
def plot_bar_horizontal(df_plot: pd.DataFrame,
                        title: str = "Contenido de GC en 16S rRNA (promedio por organismo)"):
    """
    Genera un gráfico de barras **horizontal** con %GC promedio por organismo.
    - Ordena de mayor a menor %GC.
    - Ajusta la altura de la figura según el número de organismos para que
      los nombres no queden desacomodados ni amontonados.
    """
    if "GC_promedio" in df_plot.columns:
        df_plot = df_plot.sort_values("GC_promedio", ascending=False)
        yvals = df_plot["GC_promedio"]
    else:
        df_plot = df_plot.sort_values("%GC", ascending=False)
        yvals = df_plot["%GC"]

    labels = df_plot["Organismo"]
    n = len(labels)

    # Altura dinámica: ~0.45 pulgadas por organismo
    height = max(4, 0.45 * n)
    fig, ax = plt.subplots(figsize=(10, height))

    ax.barh(labels, yvals)
    ax.set_xlabel("% GC")
    ax.set_ylabel("Organismo")
    ax.invert_yaxis()  # Organismo con mayor %GC arriba

    for i, v in enumerate(yvals):
        try:
            ax.text(float(v) + 0.3, i, f"{float(v):.1f}%", va="center", fontsize=8)
        except Exception:
            pass

    ax.set_title(title)
    ax.grid(axis="x", linestyle="--", alpha=0.4)

    plt.tight_layout()
    return fig

# -------------------------------------------------------------------
# CONCLUSIÓN GENERAL (REGLAS)
# -------------------------------------------------------------------
def conclusion_general(df_sum: pd.DataFrame) -> str:
    """
    Conclusión breve (reglas simples):
    - Mayor %GC.
    - Menor %GC.
    - Rango y promedio.
    """
    if df_sum.empty or df_sum["GC_promedio"].isna().all():
        return "Sube datos válidos para generar una conclusión."

    max_row = df_sum.loc[df_sum["GC_promedio"].idxmax()]
    min_row = df_sum.loc[df_sum["GC_promedio"].idxmin()]
    delta = float(max_row["GC_promedio"]) - float(min_row["GC_promedio"])
    mean_gc = float(df_sum["GC_promedio"].mean())

    texto = (
        f"Mayor %GC: {max_row['Organismo']} ({max_row['GC_promedio']:.2f}%).  \n"
        f"Menor %GC: {min_row['Organismo']} ({min_row['GC_promedio']:.2f}%).  \n"
        f"Rango: ~{delta:.2f} puntos porcentuales (promedio general: {mean_gc:.2f}%).\n\n"
        "En términos biológicos, un %GC elevado suele asociarse con mayor estabilidad estructural del ADN/ARN "
        "y, en algunos grupos, con tolerancia frente a ambientes exigentes. Por el contrario, un %GC reducido "
        "es frecuente en bacterias con genomas pequeños y fuerte dependencia del huésped."
    )
    return texto

# -------------------------------------------------------------------
# INTERFAZ DE USUARIO: PESTAÑAS
# -------------------------------------------------------------------
tab_carga, tab_tabla, tab_grafico, tab_conclusion = st.tabs(
    ["📥 Carga", "📊 Tabla", "📈 Gráfico", "🧠 Conclusión"]
)

# --------------------------- PESTAÑA 1: CARGA -----------------------
with tab_carga:
    st.subheader("1) Cargar archivos FASTA")
    st.caption("Formatos permitidos: .fasta, .fa, .fna (uno o varios archivos por organismo).")

    files = st.file_uploader(
        "Arrastra o selecciona uno o varios archivos FASTA",
        type=["fasta", "fa", "fna"],
        accept_multiple_files=True,
    )

    # Controles de calidad
    c1, c2, c3 = st.columns(3)
    with c1:
        len_min = st.number_input(
            "Longitud mínima 16S (bp)",
            min_value=200, max_value=3000, value=1200, step=10,
            help="Valor típico para fragmentos casi completos de 16S rRNA."
        )
    with c2:
        len_max = st.number_input(
            "Longitud máxima 16S (bp)",
            min_value=200, max_value=3000, value=1700, step=10,
            help="Límite para descartar secuencias demasiado largas o sospechosas."
        )
    with c3:
        amb_threshold = st.number_input(
            "%N máximo permitido",
            min_value=0.0, max_value=100.0, value=5.0, step=0.1,
            help="Pon 0.0 si quieres marcar cualquier N>0 como problema."
        )

    st.markdown(
        f"**Organismos base (9):** {', '.join(sorted(BASE_9))}  \n"
        f"**Adicionales sugeridos (5):** {', '.join(sorted(EXTRA_5))}"
    )

    if files:
        df_entries = build_df_from_uploads(files, len_min=len_min, len_max=len_max, amb_threshold=amb_threshold)
        df_summary = summarize_by_organism(df_entries)

        c1, c2, c3, c4 = st.columns(4)
        with c1:
            st.metric("Entradas FASTA", len(df_entries))
        with c2:
            st.metric("Organismos únicos", df_summary["Organismo"].nunique())
        with c3:
            if not df_summary.empty:
                st.metric("Promedio %GC", f"{df_summary['GC_promedio'].mean():.2f}%")
            else:
                st.metric("Promedio %GC", "N/A")
        with c4:
            if not df_summary.empty:
                rango_gc = df_summary["GC_promedio"].max() - df_summary["GC_promedio"].min()
                st.metric("Rango %GC", f"{rango_gc:.2f} pp")
            else:
                st.metric("Rango %GC", "N/A")

        bad_len = df_entries[df_entries["Flag_Longitud"]]
        bad_amb = df_entries[df_entries["Flag_Ambiguedad"]]
        st.info(
            f"Entradas totales: {len(df_entries)}  |  "
            f"Fuera de longitud: {len(bad_len)}  |  "
            f"Alta ambigüedad (≥ {amb_threshold}% N): {len(bad_amb)}"
        )

        if len(bad_len) or len(bad_amb):
            st.warning(
                "Existen entradas con problemas de longitud y/o ambigüedad. "
                "Puedes revisarlas y editar flags en la pestaña **Tabla**."
            )
        else:
            st.success("Todas las entradas pasan los criterios de longitud y ambigüedad actuales.")

        orgs_presentes = set(df_summary["Organismo"])
        st.info(
            f"Cobertura — Base: {len(BASE_9 & orgs_presentes)}/9 · "
            f"Adicionales: {len(EXTRA_5 & orgs_presentes)}/5"
        )

        st.session_state["df_entries"] = df_entries
        st.session_state["df_summary"] = df_summary

        st.success("Datos cargados y procesados. Revisa **Tabla**, **Gráfico** y **Conclusión**.")
    else:
        st.info("Aún no has subido archivos. Utiliza el control de arriba para cargar tus FASTA.")

# --------------------------- PESTAÑA 2: TABLA -----------------------
with tab_tabla:
    st.subheader("2) Tablas")

    df_entries = st.session_state.get("df_entries", pd.DataFrame())
    df_summary = st.session_state.get("df_summary", pd.DataFrame())

    if df_entries.empty:
        st.warning("Primero sube archivos en la pestaña **Carga**.")
    else:
        st.markdown("**Tabla por entrada (EDITABLE)**")

        for col in ["Flag_Longitud", "Flag_Ambiguedad"]:
            if col in df_entries.columns:
                df_entries[col] = df_entries[col].astype(bool)

        edited_entries = st.data_editor(
            df_entries,
            use_container_width=True,
            height=360,
            num_rows="dynamic",
            column_config={
                "%GC": st.column_config.NumberColumn(
                    "%GC", help="(G+C)/(A+C+G+T) * 100", step=0.001, format="%.3f"
                ),
                "%N": st.column_config.NumberColumn(
                    "%N", help="N/(A+C+G+T+N) * 100", step=0.001, format="%.3f"
                ),
                "Longitud_bases_ACGT": st.column_config.NumberColumn("Longitud ACGT"),
                "Flag_Longitud": st.column_config.CheckboxColumn(
                    "Flag_Longitud", help="Marca si la longitud NO está en el rango esperado para 16S."
                ),
                "Flag_Ambiguedad": st.column_config.CheckboxColumn(
                    "Flag_Ambiguedad", help="Marca si %N supera el umbral de ambigüedad."
                ),
            },
            key="editor_entries",
        )

        st.download_button(
            "⬇️ Descargar CSV (entradas — con ediciones)",
            edited_entries.to_csv(index=False).encode("utf-8"),
            "gc_entries.csv",
            "text/csv",
        )

        st.markdown("---")
        st.markdown("**Tabla por organismo (resumen)**")

        st.dataframe(df_summary, use_container_width=True, height=320)

        st.download_button(
            "⬇️ Descargar CSV (resumen por organismo)",
            df_summary.to_csv(index=False).encode("utf-8"),
            "gc_summary.csv",
            "text/csv",
        )

# --------------------------- PESTAÑA 3: GRÁFICO ---------------------
with tab_grafico:
    st.subheader("3) Gráfico comparativo de %GC (barras horizontales)")
    df_summary = st.session_state.get("df_summary", pd.DataFrame())

    if df_summary.empty:
        st.warning("Primero sube archivos en la pestaña **Carga** para generar el gráfico.")
    else:
        opts = list(df_summary["Organismo"])
        subset = st.multiselect(
            "Selecciona organismos a graficar (si lo dejas vacío, se grafican todos):",
            opts,
            default=opts,
        )
        df_plot = df_summary[df_summary["Organismo"].isin(subset)] if subset else df_summary

        fig = plot_bar_horizontal(
            df_plot,
            title="Contenido de GC en 16S rRNA (promedio por organismo)",
        )
        st.pyplot(fig, use_container_width=True)

        img_bytes = io.BytesIO()
        fig.savefig(img_bytes, format="png", dpi=150, bbox_inches="tight")
        st.download_button(
            "⬇️ Descargar gráfico (PNG)",
            img_bytes.getvalue(),
            "gc_barplot_horizontal.png",
            "image/png",
        )

# ------------------------ PESTAÑA 4: CONCLUSIÓN ---------------------
with tab_conclusion:
    st.subheader("4) Conclusiones automáticas")
    df_summary = st.session_state.get("df_summary", pd.DataFrame())

    if df_summary.empty:
        st.info("Sube archivos en la pestaña **Carga** para generar las conclusiones.")
    else:
        # Conclusión general breve (reglas)
        st.markdown("### 🧪 Conclusión general (reglas simples)")
        texto_general = conclusion_general(df_summary)
        st.markdown(texto_general)
        st.download_button(
            "⬇️ Descargar conclusión general (.txt)",
            texto_general.encode("utf-8"),
            "conclusion_gc_general.txt",
            "text/plain",
        )

        st.markdown("---")
        st.markdown("### 🤖 Conclusión generada por IA (OpenAI — gpt-4o-mini)")

        if st.button("Generar conclusión con IA"):
            texto_ia = generar_conclusion_openai(df_summary)
            st.session_state["conclusion_openai"] = texto_ia

        if "conclusion_openai" in st.session_state:
            st.markdown(st.session_state["conclusion_openai"])
            st.download_button(
                "⬇️ Descargar conclusión IA (.txt)",
                st.session_state["conclusion_openai"].encode("utf-8"),
                "conclusion_gc_openai.txt",
                "text/plain",
                key="download_conclusion_openai",
            )
        else:
            st.info(
                "Pulsa el botón para que el agente de IA (OpenAI) genere una conclusión extensa "
                "a partir del resumen por organismo."
            )
