# plot_descriptions.py

plot_details = {
    "1_Plot_meta": {
        "heading": "Meta Plot Overview",
        "description": "This plot provides a comprehensive overview of the meta-analysis results, highlighting key trends and patterns observed across all analyzed datasets."
    },
    "2_feature_plot2": {
        "heading": "Feature Analysis Visualization",
        "description": "Detailed visualization of key features identified in the NMF analysis, showcasing their distribution and significance across different metaprograms."
    },
    "3_split_umap_by_metaprogram": {
        "heading": "UMAP Split by Metaprogram",
        "description": "Uniform Manifold Approximation and Projection (UMAP) plot split by individual metaprograms, allowing for visual comparison of gene expression patterns across different cellular states."
    },
    "4_MP_correlation_heatmap_with_clustering": {
        "heading": "Metaprogram Correlation Heatmap",
        "description": "Heatmap displaying the correlation between different metaprograms, with hierarchical clustering to reveal relationships and similarities among cellular states."
    },
    "5_MP_Analysis_Revised_Final_treated_Donut_BArplot_Sperman": {
        "heading": "Metaprogram Analysis Summary",
        "description": "Combined visualization of metaprogram analysis results, including donut charts, bar plots, and Spearman correlation, providing a multi-faceted view of the data."
    },
    "6_pathway_enrichment_all_MPs_Hallmark": {
        "heading": "Hallmark Pathway Enrichment",
        "description": "Enrichment analysis of Hallmark gene sets across all metaprograms, highlighting key biological processes and pathways associated with each cellular state."
    },
    "7_pathway_enrichment_all_MPs_Go_BP": {
        "heading": "GO Biological Process Enrichment",
        "description": "Gene Ontology (GO) Biological Process enrichment analysis for all metaprograms, revealing functional characteristics of each identified cellular state."
    },
    "8_pathway_enrichment_all_MPs_KEGG": {
        "heading": "KEGG Pathway Enrichment",
        "description": "Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway enrichment analysis, providing insights into the molecular interaction networks associated with each metaprogram."
    },
    "9_pathway_enrichment_all_MPs_REACTOME": {
        "heading": "Reactome Pathway Enrichment",
        "description": "Enrichment analysis based on the Reactome pathway database, offering a detailed view of molecular pathways and reactions linked to each metaprogram."
    },
    "10_correlation_heatmap_enhanced": {
        "heading": "Enhanced Correlation Heatmap",
        "description": "An advanced correlation heatmap with improved visualization, showcasing the intricate relationships between various components of the analysis."
    },
    "11_correlation_corrplot": {
        "heading": "Correlation Plot",
        "description": "A correlation plot displaying the strength and direction of relationships between different variables in the study, aiding in the identification of key associations."
    },
    "12_hypergeometric_test_plot": {
        "heading": "Hypergeometric Test Results",
        "description": "Visualization of hypergeometric test results, illustrating the statistical significance of gene set overlaps and enrichments across metaprograms."
    }
}

def get_plot_details(plot_title):
    """
    Retrieve the heading and description for a given plot title.
    If not found, return default values.
    """
    details = plot_details.get(plot_title, {})
    return {
        "heading": details.get("heading", f"Plot: {plot_title}"),
        "description": details.get("description", "Details for this plot are not available.")
    }