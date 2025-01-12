API
=================


Inference
------------
.. list-table:: 
   :header-rows: 1

   * - **Main Function**
     - **Description**
   * - :doc:`get_fate_map <get_fate_map>`
     - Function for reconstructing time-scaled cell fate map.
   * - :doc:`GRNInference <GRNInference>`
     - Function for inferring dynamic gene regulatory networks along the cell fate map.


Analysis
------------
.. list-table:: 
   :header-rows: 1

   * - **Main Function**
     - **Description**
   * - :doc:`get_regulators <get_regulators>`
     - Function for retrieving the names and regulatory strengths of genes regulating a specified target gene.
   * - :doc:`get_target_genes <get_target_genes>`
     - Function for obtaining the target gene names and regulatory strengths regulated by a specified regulator gene.
   * - :doc:`find_diff_genes <find_diff_genes>`
     - Function for identifying key genes that drive the differentiation of specific progenitor cell types into downstream states.
   * - :doc:`find_fate_bias_genes <find_fate_bias_genes>`
     - Function for identifying key genes that drive the fate bias of progeny cells derived from a specific ancestral cell type.
   * - :doc:`run_fuzzy_C_means_clustering <run_fuzzy_C_means_clustering>`
     - Function for clustering regulatory interactions that coexist in similar cell types.
   * - :doc:`map_edge_clusters_to_nodes <map_edge_clusters_to_nodes>`
     - Function for identifying the specificity and constitutiveness of various clusters of regulatory interactions.




Plotting
----------
.. list-table:: 
   :header-rows: 1

   * - **Main Function**
     - **Description**
   * - :doc:`plot_dynamic_edges_number <plot_dynamic_edges_number>`
     - Function for plotting the dynamic changes in the number of regulatory interactions in the gene regulatory networks along the fate map.
   * - :doc:`plot_dynamic_target_gene <plot_dynamic_target_gene>`
     - Function for plotting the dynamic changes in the number of regulatory interactions affecting a specific target gene in the gene regulatory networks along the fate map.
   * - :doc:`plot_dynamic_regulator_number <plot_dynamic_regulator_number>`
     - Function for plotting the dynamic changes in the number of target genes regulated by a specific regulatory gene in the gene regulatory networks along the fate map.
   * - :doc:`plot_dynamic_regulator_number_heatmap <plot_dynamic_regulator_number_heatmap>`
     - Function for plotting a heatmap of the dynamic changes in the number of regulatory interactions affecting each target gene  in the gene regulatory networks along the fate map.
   * - :doc:`plot_dynamic_regulatory_network <plot_dynamic_regulatory_network>`
     - Function for plotting the dynamic changes in the regulatory network of a specific regulatory gene along the fate map.
   * - :doc:`plot_dynamic_regulator_activity <plot_dynamic_regulator_activity>`
     - Function for plotting a bar chart of the activity (number of regulated target genes) of a specific regulatory gene across different lineages and cell clusters.
   * - :doc:`plot_dynamic_regulatory_strength <plot_dynamic_regulatory_strength>`
     - Function for plotting the dynamic change in regulatory strength of each regulatory gene along a specific path of the fate map for a specific target gene.
   * - :doc:`plot_diff_genes <plot_diff_genes>`
     - Function for plotting the results of identifying key genes that drive the differentiation of specific progenitor cell types into downstream states.
   * - :doc:`plot_fate_bias_genes <plot_fate_bias_genes>`
     - Function for plotting the results of identifying key genes that drive the fate bias of progeny cells derived from a specific ancestral cell type.
   * - :doc:`plot_edge_cluster_weight <plot_edge_cluster_weight>`
     - Function for plotting the membership functions of the fuzzy clustering results of regulatory interactions, visualizing the degree of association between each regulatory interaction and its respective cluster.
   * - :doc:`plot_edge_cluster_to_nodes <plot_edge_cluster_to_nodes>`
     - Function for plotting the distribution of regulatory interactions within specific or constitutive regulatory clusters across different lineages or cell types.

