# @author Yisong Zhen
# @since  2018-09-30
# @references
#     1: Pabinger S, Dander A, Fischer M, Snajder R, Sperk M, Efremova M, Krabichler B,
#     Speicher MR, Zschocke J, Trajanoski Z. A survey of tools for variant analysis of 
#     next-generation genome sequencing data. Brief Bioinform. 2014 Mar;15(2):256-78.
#---

pkgs           <- c('tidyverse', 'DiagrammeR', 'DiagrammeRsvg')
load.packages  <- lapply(pkgs, require, character.only = TRUE)

graph <-
    create_graph() %>%
    add_n_nodes(9) %>%
    add_edges_w_string(
      '1->2 2->3 3->4 4->5
       5->6 6->7 7->8 8->9' , 'black') %>%
    set_node_attrs('label',c( 'experimental design',
                              'NGS platform illumina etc.',
                              'Quality Assesment (QC)',
                              'Read Alignment',
                              'variant Identification',
                              'Annotation',
                              'Visulization',
                              'prioritization/Fitlering',
                              'wet validation'
                             )) %>%
    select_nodes_by_id(nodes = c(1:9)) %>%
    set_node_attrs_ws('fontsize',10) %>%
    set_node_attrs_ws('fontname', 'Microsoft YaHei') %>%
    set_node_attrs_ws('fontcolor', 'black') %>%
    set_node_attrs_ws('shape','square') %>%
    set_graph_name('Roadmap')

render_graph(graph)