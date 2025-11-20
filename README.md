# Proyecto-2 - Intro_Pyhton
Este repositorio esta creado para tener un trabajo colaborativo entre los compañeros del grupo y poder sacar el proyecto adelante.

Integrantes: Roselind Vargas, Daniela Vargas y Josue Bustos. 

🧬 GC% en 16S rRNA — Dashboard Interactivo (Proyecto 2)

Dashboard interactivo desarrollado en **Streamlit** para analizar el **contenido de GC (%)** en secuencias de **16S rRNA**, permitiendo cargar archivos FASTA, comparar organismos bacterianos, generar visualizaciones, editar tablas y obtener una **conclusión automática** sobre cuál organismo presenta el mayor y el menor contenido de GC.

---

🎯 Objetivos principales

- Calcular automáticamente **%GC y %N** en secuencias 16S rRNA.
- Comparar múltiples organismos mediante tablas y gráficos.
- Responder la pregunta clave:  
  **¿Cuál organismo tiene mayor %GC y cuál menor?**
- Agregar automáticamente **5 microorganismos adicionales** al conjunto de prueba.
- Identificar fragilidades del código y proponer mitigaciones.
- Generar documentación y guion para presentación académica.
- Permitir exportar resultados en **CSV, PNG y TXT**.

---

🧪 Organismos incluidos en el proyecto

Base: 
Mycobacterium tuberculosis · Corynebacterium diphtheriae · Bacillus subtilis · Escherichia coli · Salmonella enterica · Pseudomonas aeruginosa · Clostridium botulinum · Borrelia burgdorferi · Mycoplasma genitalium

Adicionales: 
Staphylococcus aureus · Vibrio cholerae · Helicobacter pylori · Streptococcus pneumoniae · Listeria monocytogenes

---

📺 Pestañas del Dashboard

| Pestaña | Funcionalidad |
|---------|---------------|
| 📥 Carga | Subir FASTA, establecer umbrales de longitud y %N 
| 📊 Tabla | Tabla editable por archivo + resumen por organismo 
| 📈 Gráfico | Barras comparativas de %GC promedio 
| 🧠 Conclusión | Párrafo automático mayor vs menor GC + interpretación biológica 

---

🚀 Instalación y ejecución

1. Instalar dependencias

py -3.12 -m pip install streamlit pandas matplotlib

