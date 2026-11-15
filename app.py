from datetime import date
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

# Page configuration
st.set_page_config(page_title="GenoMedix", page_icon="🏥", layout="wide")

st.title("🏥 GenoMedix")
st.markdown(
    "**Clinical & Transcriptomic Informatics Suite** — Differential gene"
    " expression analysis, regional malaria surveillance, and automated"
    " diagnostic reporting."
)

# Navigation Tabs for the 3 Core Modules
tab1, tab2, tab3 = st.tabs([
    "🧬 Gene Expression Analyzer",
    "🦟 Malaria Surveillance Tracker",
    "📋 Patient Report Generator",
])

# ==========================================
# TAB 1: GENE EXPRESSION ANALYZER
# ==========================================
with tab1:
  st.subheader("Differential Gene Expression Analyzer")
  st.write(
      "Upload an RNA-seq differential expression dataset (`.csv`) or use the"
      " built-in sample dataset to identify significantly upregulated and"
      " downregulated genes."
  )

  # Option to upload CSV or use sample data
  data_source = st.radio(
      "Select Data Source:",
      ["Use Built-in Sample Dataset", "Upload Custom CSV"],
      horizontal=True,
  )

  df_expr = None

  if data_source == "Upload Custom CSV":
    uploaded_file = st.file_uploader(
        "Upload CSV (Requires columns: 'Gene', 'Log2FoldChange', 'PValue')",
        type=["csv"],
    )
    if uploaded_file is not None:
      df_expr = pd.read_csv(uploaded_file)
  else:
    # Generate realistic synthetic RNA-seq differential expression data
    np.random.seed(42)
    genes = [f"GENE_{i:03d}" for i in range(1, 151)]
    # Add some well-known clinical biomarker names
    genes[:8] = [
        "HBB",
        "TP53",
        "BRCA1",
        "IL6",
        "TNF",
        "VEGFA",
        "EGFR",
        "APOE",
    ]
    log2fc = np.random.normal(0, 1.8, 150)
    pvals = np.random.uniform(0.00001, 0.5, 150)
    # Force top biomarkers to be statistically significant
    log2fc[:4] = [3.4, -2.9, 2.5, -3.1]
    pvals[:4] = [0.0001, 0.0003, 0.001, 0.0002]

    df_expr = pd.DataFrame(
        {"Gene": genes, "Log2FoldChange": log2fc, "PValue": pvals}
    )

  if df_expr is not None:
    required_cols = {"Gene", "Log2FoldChange", "PValue"}
    if not required_cols.issubset(df_expr.columns):
      st.error(
          "CSV must contain the exact columns: `Gene`, `Log2FoldChange`, and"
          " `PValue`."
      )
    else:
      # Interactive Threshold Sliders
      col_t1, col_t2 = st.columns(2)
      fc_cutoff = col_t1.slider("Log2 Fold Change Cutoff", 0.5, 4.0, 1.5, 0.1)
      pval_cutoff = col_t2.number_input(
          "P-Value Significance Threshold",
          min_value=0.0001,
          max_value=0.1,
          value=0.05,
          format="%.4f",
      )

      # Classify genes
      df_expr["-Log10PValue"] = -np.log10(df_expr["PValue"])

      def classify_gene(row):
        if (
            row["PValue"] < pval_cutoff
            and row["Log2FoldChange"] >= fc_cutoff
        ):
          return "Upregulated"
        elif (
            row["PValue"] < pval_cutoff
            and row["Log2FoldChange"] <= -fc_cutoff
        ):
          return "Downregulated"
        return "Not Significant"

      df_expr["Regulation"] = df_expr.apply(classify_gene, axis=1)

      # Summary Metrics
      up_count = (df_expr["Regulation"] == "Upregulated").sum()
      down_count = (df_expr["Regulation"] == "Downregulated").sum()

      m1, m2, m3 = st.columns(3)
      m1.metric("Total Genes Analyzed", len(df_expr))
      m2.metric("Upregulated Biomarkers", up_count)
      m3.metric("Downregulated Biomarkers", down_count)

      # Interactive Volcano Plot
      fig = px.scatter(
          df_expr,
          x="Log2FoldChange",
          y="-Log10PValue",
          color="Regulation",
          hover_name="Gene",
          color_discrete_map={
              "Upregulated": "#ef4444",
              "Downregulated": "#3b82f6",
              "Not Significant": "#9ca3af",
          },
          title="Volcano Plot: Differential Gene Expression",
      )
      fig.add_hline(
          y=-np.log10(pval_cutoff), line_dash="dash", line_color="gray"
      )
      fig.add_vline(x=fc_cutoff, line_dash="dash", line_color="gray")
      fig.add_vline(x=-fc_cutoff, line_dash="dash", line_color="gray")
      st.plotly_chart(fig, use_container_width=True)

      # Filterable Data Table
      st.write("### Significant Targets")
      sig_df = df_expr[df_expr["Regulation"] != "Not Significant"].sort_values(
          by="PValue"
      )
      st.dataframe(
          sig_df[["Gene", "Log2FoldChange", "PValue", "Regulation"]],
          use_container_width=True,
      )

# ==========================================
# TAB 2: MALARIA SURVEILLANCE TRACKER
# ==========================================
with tab2:
  st.subheader("🦟 Regional Malaria Surveillance & Prevalence Tracker")
  st.write(
      "Monitors *Plasmodium falciparum* incidence rates, rapid diagnostic test"
      " (RDT) positivity, and antimalarial resistance markers across sentinel"
      " sites."
  )

  # Sample epidemiological dataset
  malaria_data = pd.DataFrame({
      "Sentinel Site": [
          "Port Harcourt Central",
          "Obio-Akpor Clinic",
          "Eleme General",
          "Bonny Island Hub",
          "Ahoada West",
          "Degema Station",
      ],
      "lat": [4.8156, 4.8920, 4.7924, 4.4342, 5.0828, 4.7495],
      "lon": [7.0498, 6.9922, 7.1213, 7.1706, 6.6498, 6.7689],
      "Cases_Reported": [420, 510, 290, 185, 340, 215],
      "RDT_Positivity_Rate (%)": [34.5, 41.2, 28.0, 19.4, 38.9, 25.1],
      "PfKelch13_Resistance (%)": [2.1, 3.4, 1.5, 0.8, 2.9, 1.2],
  })

  # Top KPI Metrics
  k1, k2, k3 = st.columns(3)
  k1.metric(
      "Total Confirmed Cases",
      f"{malaria_data['Cases_Reported'].sum():,}",
  )
  k2.metric(
      "Mean RDT Positivity",
      f"{malaria_data['RDT_Positivity_Rate (%)'].mean():.1f}%",
  )
  k3.metric(
      "Avg Artemisinin Resistance Marker",
      f"{malaria_data['PfKelch13_Resistance (%)'].mean():.1f}%",
  )

  col_map, col_chart = st.columns([1.2, 1])

  with col_map:
    st.write("### Geospatial Distribution")
    st.map(malaria_data, latitude="lat", longitude="lon", size="Cases_Reported")

  with col_chart:
    st.write("### Positivity Rate by Sentinel Site")
    bar_fig = px.bar(
        malaria_data.sort_values("RDT_Positivity_Rate (%)", ascending=True),
        x="RDT_Positivity_Rate (%)",
        y="Sentinel Site",
        orientation="h",
        color="RDT_Positivity_Rate (%)",
        color_continuous_scale="Reds",
    )
    st.plotly_chart(bar_fig, use_container_width=True)

  st.write("### Sentinel Site Surveillance Log")
  st.dataframe(
      malaria_data.drop(columns=["lat", "lon"]), use_container_width=True
  )

# ==========================================
# TAB 3: PATIENT REPORT GENERATOR
# ==========================================
with tab3:
  st.subheader("📋 Automated Clinical Diagnostic Report Generator")
  st.write(
      "Enter patient laboratory findings to generate a standardized diagnostic"
      " report ready for export."
  )

  with st.form("patient_report_form"):
    c1, c2, c3 = st.columns(3)
    patient_name = c1.text_input("Patient Full Name", value="Chinedu Okafor")
    patient_id = c2.text_input("Patient / Lab ID", value="GMX-2026-084")
    age_sex = c3.text_input("Age / Sex", value="29 / Male")

    c4, c5, c6 = st.columns(3)
    genotype = c4.selectbox(
        "Hemoglobin Genotype (HBB)", ["HbAA", "HbAS", "HbSS", "HbAC", "HbSC"]
    )
    malaria_rdt = c5.selectbox(
        "Malaria RDT / Microscopy",
        [
            "Negative",
            "Positive (P. falciparum +)",
            "Positive (P. falciparum +++)",
        ],
    )
    hb_level = c6.number_input(
        "Hemoglobin Level (g/dL)",
        min_value=3.0,
        max_value=20.0,
        value=13.5,
        step=0.1,
    )

    clinician_notes = st.text_area(
        "Clinical & Molecular Interpretation Notes",
        value=(
            "Patient screened for hemoglobinopathy and malaria parasitemia."
            " Recommend routine follow-up."
        ),
    )

    generate_btn = st.form_submit_button("Generate Clinical Report")

  if generate_btn:
    # Automated clinical flagging logic
    anemia_status = (
        "Normal" if hb_level >= 12.0 else "Low (Anemia Indicated)"
    )

    report_text = f"""============================================================
GENOMEDIX CLINICAL & MOLECULAR DIAGNOSTIC REPORT
============================================================
Report Date:       {date.today().strftime('%B %d, %Y')}
Patient ID:        {patient_id}
Patient Name:      {patient_name}
Age / Sex:         {age_sex}
------------------------------------------------------------
LABORATORY & MOLECULAR FINDINGS
------------------------------------------------------------
1. Hemoglobin Genotype (HBB):   {genotype}
2. Malaria Parasite Screen:     {malaria_rdt}
3. Hemoglobin Concentration:    {hb_level:.1f} g/dL [{anemia_status}]
------------------------------------------------------------
CLINICAL INTERPRETATION
------------------------------------------------------------
{clinician_notes}

============================================================
Generated via GenoMedix Clinical Informatics Suite
============================================================"""

    st.success("Diagnostic report compiled successfully.")
    st.code(report_text, language="text")

    # Direct Download Button
    st.download_button(
        label="📥 Download Report (.txt)",
        data=report_text,
        file_name=f"{patient_id}_Diagnostic_Report.txt",
        mime="text/plain",
  )
    
