# Introduction

This page provides notes and instructions on how to explore, visualize and extract drug-similarity information in DrugSimDB.

Current version of DrugSimDB includes chemical, structure, pathway, and functional GO terms similarities of <b>10,317</b> small-molecule drugs, yielding <b>238,635</b> drug-pairs with significant aggregated score.

Overall, this interface is user-friendly and responsive.


<img style='display: table; border-radius: 5px; border: 1px solid #293954; '
              src="documentation/Figure 1.png" width="400">

# Main Search Interface
The main interface for exploring the DrugSimDB database looks like bellow:

<img style='display: table; border-radius: 5px; border: 1px solid #293954; '
              src="documentation/Figure 2.png" width="400">

## Search a Drug
The following image shows the search panel, where a drug name is queried in the DrugSimDB database to retrieve its similarity information with other drugs. Please note DrugSimDB only includes similarity information about statistically significant drug-pairs (<b>total: 238,635</b>). 

<!-- ![](../www/tutorial images/search_pan.png) -->
<img style='display: table; border-radius: 5px; border: 1px solid #293954; '
              src="documentation/search_pan.png" width="300">

It lists all the available drugs that have similarity information in the DrugSimDB, including a search bar and an option for choosing the combined similarity score to display with, i.e. aggregated score, p-value or adjusted p-value. The search bar supports auto-complete feature (see below). Therefore, when user enters a partial keyword as the drug name, for example, "<b>biva</b>", the list automatically filters to matching drug names "<b>Bivalirudin</b>" and "<b>Tetrahydrocannabivarin</b>". 

<img style='display: table; border-radius: 5px; border: 1px solid #293954; '
              src="documentation/search_pan_search_demo.png" width="300">

Each Drug name has a "<b>View</b>" button next to it, and when clicked, the result is displayed in the right panel that includes all the similarity information about the query drug (see next sub-section).

### Query Result (Tabular-View)
When a the <b>View</b> button is pressed for a drug, all other drugs that have statistically significant combined-similarity will be shown in a table along with other individual similarity scores, e.g. structure similarity, chemical similarity, target pathway similarity, GO terms similarity. The table also shows the statuses of the incident drugs annotated in [Drugbank](https://www.drugbank.ca). Note, depending on the option chosen (see above) before the query being made, the final column in the table displays the information about the combined knowledge of all the individual similarity scores, i.e. Combined score, p-value, or adjusted p-value. 

<img style='display: table; border-radius: 5px; border: 1px solid #293954; '
              src="documentation/search_tableView.png" width="600">

This table view is exportable in multiple format, e.g. CSV, Excel, PDF. It can also be copied in the ClipBoard and be printed using system default printer driver. There is a global search bar at the top of the view, which allows user to filter any the whole table based on a query word that works for all the columns. Morever, the table can be sorted based on individual column and has pagination support.

### Query Result (Network-view)
One of the most important feature of this interface is the network view of drug-pairs based on the similarity. Including the query drug, a subset of all the other drugs that are reported as being significantly similar (see the Table-View) in the DrugSimDB, has been used to extract the induced sub-network from the DrugSimDB network (undirected, formed by all the drug-pairs). The network-view of that induced sub-network looks like the following:


<img style='display: table; border-radius: 5px; border: 1px solid #293954; '
              src="documentation/search_NetworkView.png" width="600">

This view has several useful feature, which are described below.

#### <b>Interactive</b>
This view is interactive. Users can drag the components of the network, i.e. nodes or edges using mouse or touch in order for better visualization. Moreover, the view-pan has navigation button to move the whole network-view left, right, up, or down. There are buttons for zooming-in, zooming-out and expanding the network to fit the pan. Usuers can also export the whole view-pan, which only supports .png format only. 

#### <b>Visualy informative</b>
IN this network-view, the nodes (Drugs) are colored based on the groups that they were annotated in the [Drugbank](https://www.drugbank.ca), e.g. <b>approved</b>; <b>approved, experimental</b>; and <b>approved, vet_approved</b>. All the edge colors are the same, but their widths varies in proportionate to the scores shown in the final column (combined score, or p-value, or adjusted p-value) of the tablur-view. 

#### <b>Query-able</b>
This network-view is query-able. 

- <b><i>Node:</i></b> On node mouse-hover event, user can see the following informaiton 
    1. Drug Name
    2. Status, i.e. Approved, investigational, etc.
    3. List of adverse Side-effects (from [SIDER](http://sideeffects.embl.de/) database)
    
    
    <img style='display: table; border-radius: 5px; border: 1px solid #293954; '
              src="documentation/node_hover.png" width="300">
    
- <b><i>Edge:</i></b> On click-event of an edge, it makes a pubmed query (using [easyPubMed R package](https://cran.r-project.org/web/packages/easyPubMed/easyPubMed.pdf)) with the <b>two incident drug names</b>, and retrieves and p a list of pubmed articles that have the co-occurance of those two drugs. This process has several steps:

    - <i>Step #1</i>: Using both drug names (incident to an edge), make a query in PubMed in all the fields, and a list of article IDs is returned, where the drug-pair is co-mentioned in corresponding articles. For example, for an edge with drug-pair, <b>Apixaban-Bivalirudin</b>, the PubMed query string forms like this: <b>Apixaban[All Fields] AND Bivalirudin[All Fields]</b>. We have used the <i>get_pubmed_ids()</i> function defined in the [easyPubMed R package](https://cran.r-project.org/web/packages/easyPubMed/easyPubMed.pdf) for this purpose.
    
    - <i>Step #2:</i>Next, with the article IDs, a set of pubmed article data is retrieved from PubMed by using  <i>fetch_pubmed_data()</i> function defined in the [easyPubMed R package](https://cran.r-project.org/web/packages/easyPubMed/easyPubMed.pdf). This step is currently restricts fetching PumMed data for top 20 records in the list. Hence, even if the number of found article IDs can be arbitrarily large, this constraint keeps the network and cpu usage within tractable scope.
    
    - <i>Step #3:</i> The <i>fetch_pubmed_data()</i> function returns data in XML format which requires further processing, i.e. converting into list () and then finally a dataframe of article info (using our in-house script), using using the <i>articles_to_list()</i> and <i>article_to_df()</i> functions defined in the [easyPubMed R package](https://cran.r-project.org/web/packages/easyPubMed/easyPubMed.pdf). This dataframe contains article information which then later displayed in modal window (see below) in the PubMed reference format. Note, since the dataframe is tabularied as a data table, user can 1) search the that table using a global search bar at the top-right cornar of that modal window, 2) navigate the list using the pagination shown at the bottom of the table.
    
        <img style='display: table; border-radius: 5px; border: 1px solid #293954; '
              src="documentation/pubmedSearch.png" width="700">


### Query drug info
For the query drug, it also provides following information:
  1. P


### Query drug info
TBA

## Download
TBA

## Statistics
TBA

## Issues/Bugs
TBA

