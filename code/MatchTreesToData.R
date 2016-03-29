##########################################################################################################
# Functions for pruning lgfam trees (Dediu 2015, https://github.com/ddediu/lgfam-newick) 
# so that the tips match the data we have, based on an ID of choice (ISO, AUTOTYP or GLOTTOLOG IDs)
# 
# Because data can be collected at various taxonomic levels (e.g. at the dialect or the language level),
# it sometimes happens that we have data on nodes rather than tips in a tree, e.g. we might have data 
# for German but our tree ends in tips that represent only the dialects German. 
# This creates a problem for all functions that expect data to match tips only (e.g. virtually all
# phylogenetic methods in `ape' or `geiger'). 
# To maximize data usage, we fill the tips that lack data with the ID of the next higher node. 
# This is justified on the assumption that one would have gathered data from the daugthers of this node 
# if one suspected strong variation. 
# In the absence of this, it is reasonable to assume that the data of the mother node are also 
# characteristic of the daughters. This idea is implemented in the function fill.tree() which takes 
# a list of trees or a multiPhylo object, based on various language identifiers, as specified by the 
# "id.type" argument
#
# The second function, prune.to.data(), prunes the trees to retain only tips with data.
#
# MAY CONTAIN BUGS. USE WITH CARE!
#
# Balthasar Bickel [2015-11-16 BB]
#
# If you use this code, please cite it!
# 
##########################################################################################################

require(dplyr)
require(ape)
require(phylobase)

fill.tree <- function(trees, ids, id.type = c('ISO', 'AUTOTYP.LID', 'GLOTTOCODE'), ...) {
	
	# tree is a multiPhylo object or a list of trees
	# ids is a vector of IDs that index the data we have
	# id.type is the choice of ID type
	
	# -- choose ID type:
	id.regexes <- list(
			ISO = '.*\\[i\\-([a-z]{3,4})\\].*', 	
	        	AUTOTYP.LID = '.*\\[a\\-(\\d{1,4})\\].*',
			GLOTTOCODE = '.*\\[g\\-([a-z]{4}\\d{4})\\].*'
			)			
	id.regex <- id.regexes[[match.arg(id.type)]]
	

	# -- search for ID codes among nodes and use them in the tips:
	lapply(trees, function(t) {
	
		# make all tips that have ID codes pure IDs:
		t$tip.label <- ifelse(  grepl(id.regex, t$tip.label, perl=T),
					paste(gsub(id.regex, "\\1", t$tip.label, perl=T)),
					paste(t$tip.label)
					)

		# get the nodes which have IDs (we do this without transforming node labels to IDs for fear of empty node labels)				
		id.nodes <- grep(id.regex, t$node.label, value=T, perl=T)

		# but keep only those that we actually also have data on:
		id.nodes <- id.nodes[gsub(id.regex, '\\1', id.nodes, perl=T) %in% ids] 
			# print(id.nodes)

		# get the tips of the ID-bearing nodes 
		for (node in id.nodes) {

			# extract all tips of the ID-bearing node without those that have data (because it can happen that the same ID occurs both as a node and a tip, in which case we keep the tip because it is more informative about branch lengths and topology)
			# extract.clade only works if there at least 2 daughters left. Instead, we use phylobase functions:
			tips <- names(suppressWarnings(descendants(phylo4(t), node=node, type='tips')))
				# print(node)
				# print(tips)		
				# clade <- extract.clade(t, node)
				# print(clade$tip.label)
				# empty.tips <- clade$tip.label[!clade$tip.label  %in% ids]
			empty.tips <- tips[!tips %in% ids]

			# now declare the first (random) tip as a data-bearing node:
			t$tip.label <- ifelse(t$tip.label %in% empty.tips, 
					paste(gsub(id.regex, '\\1', node, perl=T)), 
					paste(t$tip.label)
					)
 			}
		# in some cases, this creates duplicates, i.e. tips that now have the same name. We pick one randomly or, for simplicity, whatever happens to be the first. But we first need to catch trees that contain only duplicates (e.g. Basque, Anson Bay or Eastern Daly). These would result in isolates, not amenable to standard tree tests:
		if (length(unique(t$tip.label))==1) {
			t <- NULL
			} else {			
		t <- drop.tip(t, which(duplicated(t$tip.label)))
			}
		return(t)
		}
		) %>% 
	.[!sapply(., is.null)] # flushing out trees that got reduced to singletons/isolates because all tips have the same ID
	}		

prune.to.data <- function(trees, ids, ...) {
	
	# tree is a multiPhylo object or a list of trees. Tips must have the same format (e.g. ISO or glottocode) as the ids
	# ids is a vector of ID which index data we have. They must be in the same format as the tip.label in the tree.
	
	lapply(trees, function(t) {
		# prune so that we only have tips for which we have data. 
		# drop.tip throws an error if the pruning leaves a tree with only one tip, so we prune only if we will get a tree with at least 2 tips and discard isolates since we can't use them for standard tree tests anyway:
		if (length(intersect(t$tip.label, ids))>1) {
			t <- drop.tip(t, setdiff(t$tip.label, ids))
			} else {t <- NULL}
	}) %>% 
	.[!sapply(., is.null)] # flushing out empty trees with no data
	}