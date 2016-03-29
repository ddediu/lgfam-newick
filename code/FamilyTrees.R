###############################################################################################################
#
# This script implements language family tres and collections thereof.
#
# Trees are implemented by the familytree class, an extension of class phylo,
# and collections of trees (a classification) by the languageclassification class.
#
# The node names are composed of the linguistic entity's name (as it appears in the classification),
# followed by a code of the form [i-ISO][w-WALS][a-AUTOTYP][g-GLOTTOCODE], 
# where any of the codes can be missing (but the symbol must be present, e.g. [i-iso][w-][a-][g-glottocode])
# or can appear multiple times separated by dash (e.g. [i-iso][w-wals1-wals2][a-autotyp][g-glottocode]).
#
# Copyright (C) 2013-2015  Dan Dediu
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#################################################################################################################


library(parallel);  # multicore processing
library(stringr);   # string processing
library(compiler);  # compiler for better speed
library(phytools);  # phylo class and read.newick()
library(phangorn);  # nnls.tree()
library(GA);        # genetic algorithms

# Print error message?
PRINT_DEBUG_MESSAGES = FALSE;

# GA parameters:
GA.MAXITER     = 10000; # maximum number of iterations to run
GA.POPSIZE     = 100;   # population size
GA.CONSTANTRUN = 100;   # the number of consecutive generations resulting in the same fitness needed to stop the search prematurely


###########################################
#
#   Language and family name processing
#
###########################################

# Replace characters by other characters in a string:
.replace.multiple.chars <- function(string, replacements=NULL)
{
  if( is.null(string) || is.na(string) || string=="" || is.null(replacements) || is.na(replacements) || replacements=="" ) return (string);
  for( i in 1:length(replacements) ) string <- gsub( names(replacements)[i], replacements[i], string, fixed=TRUE );
  return (string);
}
.replace.multiple.chars <- cmpfun(.replace.multiple.chars);

# Replace troublesome characters in language names:
normalize.language.name <- function( string )
{
  return (.replace.multiple.chars(string, 
                                  c(# ASCII characters not accepted by the Newick format:
                                    ","=".", "'"="`", "\\'"="`", "’"="`", "("="{", ")"="}", "\t"=" ", ":"="|", ";"="|",
                                    # Non-ASCII and combinations that must be converted to closesnt ASCII:
                                    "á"="a", "à"="a", "ã"="a", "ä"="a", "â"="a", "Ä"="A", "Á"="A", "À"="A",
                                    "é"="e", "è"="e", "ë"="e", "ê"="e", "É"="E", 
                                    "í"="i", "ì"="i", "î"="i", "ï"="i", 
                                    "ó"="o", "ò"="o", "ö"="o", "õ"="o", "ô"="o", "Ö"="O",
                                    "ü"="u", "ù"="u", "ú"="u",
                                    "ç"="c",
                                    "ñ"="n", "Ñ"="N"
                                    )
                                  ));
}
normalize.language.name <- cmpfun(normalize.language.name);
# Check a set of roots for weird characters:
# cat("",file="names.txt",append=FALSE);for(i in 1:length(roots)){ cat(paste(unique(vapply(c(extract.languages(roots[[i]]), extract.internal.nodes(roots[[i]])),function(s) extract.name.and.codes(s)$name, character(1))), collapse="\n"),file="names.txt",append=TRUE); cat("\n",file="names.txt",append=TRUE)}

# Make sure strings are in the correct format for the Newick/Nexus formats:
nexus.string <- function( string, force.to.nexus=FALSE )
{
  if( force.to.nexus )
  {
    return (.replace.multiple.chars(string, c(" "="_", "-"="_", ","="_", "("="{", ")"="}", "/"="", "'"="`", ":"="_", ";"="_", "\""="", "+"="_and_", "."="")));
  } else
  {
    return (string);
  }
}
nexus.string <- cmpfun(nexus.string);

# For a given qualified name (possibly quoted and with codes [i-][w-][a-][g-]), extract the proper name and the codes:
# The glottolog.convention concerns the language codes: it gives [glottocode]{[iso]} where {}=optional:
extract.name.and.codes <- function( full.name, quotes="'", glottolog.convention=FALSE, language.name.can.be.empty=FALSE )
{
  if( !is.na(quotes) )
  {
    full.name <- substr( full.name, 2, nchar(full.name)-1 );
  }
  
  if( is.null(full.name) || is.na(full.name) || full.name=="" )
  {
    if( PRINT_DEBUG_MESSAGES ) cat( "Empty name!\n" );
    return ( list( "name"=full.name, "iso"="", "wals"="", "autotyp"="", "glottolog"="" ) );
  }
  
  # Extract the actual name:
  name <- str_trim( strsplit( full.name, "[", fixed=TRUE )[[1]][1] );
  if( is.null(name) || name=="" )
  {
    if( !language.name.can.be.empty )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "Malformed qualified name: the actual name is not specified!\n" );
      return (NULL);
    } else
    {
      name <- "";
    }
  }
  
  isocode <- walscode <- autotypcode <- glottocode <- "";
  if( !glottolog.convention )
  {
    # Extract the iso code(s):
    isocode <- strsplit( full.name, "[i-", fixed=TRUE )[[1]];
    if( length(isocode) < 2 )
    {
      # Iso code not given:
      isocode <- NULL;
    } else
    {
      isocode <- strsplit( isocode[2], "]", fixed=TRUE )[[1]]
      isocode <- isocode[1];
      if( length(grep("-",isocode,fixed=TRUE)) > 0 ) isocode <- strsplit( isocode, "-", fixed=TRUE )[[1]];
    }
    
    # Extract the wals code:
    walscode <- strsplit( full.name, "[w-", fixed=TRUE )[[1]];
    if( length(walscode) < 2 )
    {
      # wals code not given:
      walscode <- NULL;
    } else
    {
      walscode <- strsplit( walscode[2], "]", fixed=TRUE )[[1]]
      walscode <- walscode[1];
      if( length(grep("-",walscode,fixed=TRUE)) > 0 ) walscode <- strsplit( walscode, "-", fixed=TRUE )[[1]];
    }
    
    # Extract the autotyp code:
    autotypcode <- strsplit( full.name, "[a-", fixed=TRUE )[[1]];
    if( length(autotypcode) < 2 )
    {
      # autotyp code not given:
      autotypcode <- NULL;
    } else
    {
      autotypcode <- strsplit( autotypcode[2], "]", fixed=TRUE )[[1]]
      autotypcode <- autotypcode[1];
      if( length(grep("-",autotypcode,fixed=TRUE)) > 0 ) autotypcode <- strsplit( autotypcode, "-", fixed=TRUE )[[1]];
    }
    
    # Extract the glottolog code:
    glottocode <- strsplit( full.name, "[g-", fixed=TRUE )[[1]];
    if( length(glottocode) < 2 )
    {
      # glottolog code not given:
      glottocode <- NULL;
    } else
    {
      glottocode <- strsplit( glottocode[2], "]", fixed=TRUE )[[1]]
      glottocode <- glottocode[1];
      if( length(grep("-",glottocode,fixed=TRUE)) > 0 ) glottocode <- strsplit( glottocode, "-", fixed=TRUE )[[1]];
    }
  } else
  {
    # Extract the glottocode and the optional iso code:
    codes <- str_trim( strsplit( full.name, "[", fixed=TRUE )[[1]][-1] );
    if( length(codes) > 0 )
    {
      # There's at least the glottocode:
      glottocode <- strsplit( codes[1], "]", fixed=TRUE )[[1]][1];
      if( length(codes) > 1 )
      {
        # There's also an iso:
        isocode <- strsplit( codes[2], "]", fixed=TRUE )[[1]][1];
      }
    }
  }
  
  return ( list( "name"=name, "iso"=isocode, "wals"=walscode, "autotyp"=autotypcode, "glottolog"=glottocode ) );
}
extract.name.and.codes <- cmpfun(extract.name.and.codes);

# Create a qualified language name given the codes (possibly quoted and with codes [i-][w-][a-][g-]):
create.name.and.codes <- function(name, iso="", wals="", autotyp="", glottocode="", quotes="'" )
{
  paste0(quotes, normalize.language.name(name), 
         " [i-", paste0(iso,collapse="-"), "][w-", paste0(wals,collapse="-"), "][a-", paste0(autotyp,collapse="-"), "][g-", paste0(glottocode,collapse="-"), "]", 
         quotes);
}
create.name.and.codes <- cmpfun(create.name.and.codes);

# Test if two nodes contain the same language info:
are.nodes.equal <- function( x, y, quotes=NA, are.nodes.parsed=FALSE )
{
  if( !are.nodes.parsed )
  {
    x <- extract.name.and.codes(x, quotes=quotes, language.name.can.be.empty=TRUE);
    y <- extract.name.and.codes(y, quotes=quotes, language.name.can.be.empty=TRUE);
  }
  
  # Helper function: any string is equal to NA or "", but if defined they must be identical:
  .helper.equal <- function( s1, s2 )
  {
    if( is.na(s1) || s1 == "" || is.na(s2) || s2 == "" )
    {
      return (TRUE); # whatever...
    } else
    {
      return (s1 == s2); # but if defined they must really be equal
    }
  }
  
  # Helper function: check the equality of two (lists) of codes, taking into account permutations and inclusions:
  .helper.equal.codes <- function( c1, c2, strictly.equal=FALSE )
  {
    if( is.null(c1) || is.na(c1) || c1 == "" || is.null(c2) || is.na(c2) || c2 == "" )
    {
      return (TRUE); # whatever...
    } else
    {
      # But if defined they must really be equal
      return ( ifelse( strictly.equal, setequal(c1,c2), length( intersect(c1,c2) ) > 0 ) ); 
    }
  }
  
  return ( .helper.equal( x$name, y$name ) && 
             .helper.equal.codes( x$iso, y$iso ) && 
             .helper.equal.codes( x$wals, y$wals ) && 
             .helper.equal.codes( x$autotyp, y$autotyp ) && 
             .helper.equal.codes( x$glottolog, y$glottolog ) );
}
are.nodes.equal <- cmpfun(are.nodes.equal);

# Convert from the glottolog to my convention filling in the extra codes as well:
.glottolog2mine <- function(name, mapping, quotes="'")
{
  if( is.null(name) || is.na(name) || name=="" ) return (name);
  # Extract the codes from the glottolog convention:
  tmp <- extract.name.and.codes(name, quotes, glottolog.convention=TRUE);
  if( is.null(tmp) ) return (name);
  # Fill in the missing codes:
  if( !is.null(mapping) )
  {
    # Try to fill in the missing codes, attempting to match the given information:
    s <- NULL;
    if( length(s) == 0 && (length(tmp$glottolog) > 1 || tmp$glottolog != "") ) s <- which(mapping$Glottolog %in% tmp$glottolog);
    if( length(s) == 0 && (length(tmp$iso) > 1       || tmp$iso != "")       ) s <- which(mapping$ISO       %in% tmp$iso      );
    if( length(s) == 0 && (length(tmp$wals) > 1      || tmp$wals != "")      ) s <- which(mapping$WALS      %in% tmp$wals     );
    if( length(s) == 0 && (length(tmp$autotyp) > 1   || tmp$autotyp != "")   ) s <- which(mapping$AUTOTYP   %in% tmp$autotyp  );
    s <- unique(s);
    if( length(s) > 0 )
    {
      # Found mappings: use them!
      if( length(tmp$iso) == 1       && tmp$iso       == "" ){ tmp$iso       <- unique(mapping$ISO[s]      ); tmp$iso       <- tmp$iso      [ tmp$iso      != "" ]; }
      if( length(tmp$wals) == 1      && tmp$wals      == "" ){ tmp$wals      <- unique(mapping$WALS[s]     ); tmp$wals      <- tmp$wals     [ tmp$wals     != "" ]; }
      if( length(tmp$autotyp) == 1   && tmp$autotyp   == "" ){ tmp$autotyp   <- unique(mapping$AUTOTYP[s]  ); tmp$autotyp   <- tmp$autotyp  [ tmp$autotyp  != "" ]; }
      if( length(tmp$glottolog) == 1 && tmp$glottolog == "" ){ tmp$glottolog <- unique(mapping$Glottolog[s]); tmp$glottolog <- tmp$glottolog[ tmp$glottolog != "" ]; }
    }
  }
  iso      <- paste( tmp$iso,      collapse="-", sep="" );
  wals     <- paste( tmp$wals,     collapse="-", sep="" );
  autotyp  <- paste( tmp$autotyp,  collapse="-", sep="" );
  glottolog <- paste( tmp$glottolog, collapse="-", sep="" );
  
  # Convert them to the extended notation:
  return (paste0( ifelse(!is.na(quotes),quotes,""), tmp$name, " ", "[i-", iso, "][w-", wals, "][a-", autotyp, "][g-", glottolog, "]", ifelse(!is.na(quotes),quotes,"") ));
}
.glottolog2mine <- cmpfun(.glottolog2mine);


###########################################################################
#
# A tree is an extension of the phylo class
#
###########################################################################

# Create a familytree object from a Newick representation (in a file or as text):
familytree <- function(text=NULL, tree.name="", file.name="")
{
  if( (is.null(text) || is.na(text) || text=="") && (file.name == "" || is.null(file.name) || is.na(file.name)) )
  {
    warning("To build a 'familytree' object a 'text' or 'file' must be given!\n");
    return (NULL);
  }
  # Read the Newick format:
  p <- read.newick(file.name, text);
  if( is.null(p) ) return (NULL);
  # Create the familytree object:
  attr(p, "tree.name") <- tree.name;
  class(p) <- append("familytree", class(p));
  return (p);
}
# tree <- familytree("(((Chicomuceltec[i-cob][w-cec][a-1167][g-chic1271]:0.777778)Huastecan[i-][w-][a-][g-huas1241]:0.111111,(((((Chol[i-ctu][w-][a-1196][g-chol1282]:0.333333,(Buena_Vista_Chontal[i-][w-][a-][g-buen1245]:0.222222,Miramar_Chontal[i-][w-][a-][g-mira1253]:0.222222,Tamulté_de_las_Sábanas_Chontal[i-][w-][a-][g-tamu1247]:0.222222)Tabasco_Chontal[i-chf][w-cmy][a-1136][g-taba1266]:0.111111)Chol-Chontal[i-cti][w-col][a-][g-chol1281]:0.111111,(Cholti[i-][w-][a-][g-chol1283]:0.333333,Chorti[i-caa][w-coi][a-1105][g-chor1273]:0.333333)Chorti-Cholti[i-][w-][a-][g-chor1272]:0.111111,Epigraphic_Mayan[i-emy][w-][a-][g-epig1241]:0.444444)Cholan[i-][w-][a-][g-chol1287]:0.111111,((Chanal_Cancuc[i-][w-][a-][g-chan1320]:0.333333,Tenango[i-][w-][a-][g-tena1239]:0.333333)Tzeltal[i-tzh][w-tza-tzt][a-1433-1434-2651-804][g-tzel1254]:0.111111,Tzotzil[i-tzo][w-][a-2652-2654-2655][g-tzot1259]:0.444444)Tzeltalan[i-tzb][w-tze][a-][g-tzel1253]:0.111111)Cholan-Tzeltalan[i-][w-][a-][g-chol1286]:0.111111,((Chuj[i-cac][w-][a-1107][g-chuj1250]:0.444444)Chujean[i-cac][w-chj][a-][g-chuj1249]:0.111111,((Popti`[i-jac][w-jak][a-460][g-popt1235]:0.333333,Q`anjob`al[i-kjb][w-kea][a-1760][g-qanj1241]:0.333333,Western_Kanjobal[i-knj][w-kwe][a-1787][g-west2635]:0.333333)Kanjobal-Jacaltec[i-][w-][a-][g-kanj1263]:0.111111,(Motozintleco[i-][w-][a-][g-moto1243]:0.333333,Tuzanteco[i-][w-][a-][g-tuza1238]:0.333333)Mocho[i-mhc][w-][a-][g-moch1257]:0.111111)Kanjobalan[i-][w-][a-][g-kanj1262]:0.111111)Kanjobalan-Chujean[i-][w-][a-][g-kanj1261]:0.111111,(((Aguacateco[i-agu][w-agu][a-855][g-agua1252]:0.333333,Ixil[i-ixl][w-][a-1665][g-ixil1251]:0.333333)Ixilan[i-ixi][w-ixi][a-][g-ixil1250]:0.111111,(Mam[i-mam][w-][a-2039-2101][g-mamm1241]:0.333333,Tektiteko[i-ttc][w-tec][a-2627][g-tekt1235]:0.333333)Mamean[i-][w-][a-][g-mame1240]:0.111111)Greater_Mamean[i-][w-][a-][g-grea1277]:0.111111,(Kekchi[i-kek][w-kek][a-1730][g-kekc1242]:0.444444,(((Kaqchikel[i-cak][w-][a-][g-kaqc1270]:0.111111,Tz`utujil[i-tzj][w-][a-388][g-tzut1248]:0.111111)Cakchiquel-Tzutujil[i-tzj][w-tzu][a-][g-cakc1244]:0.111111,(Achi[i-acr][w-][a-826][g-achi1256]:0.111111,K`iche`[i-quc][w-qch][a-337][g-kich1262]:0.111111)Quiche-Achi[i-acc][w-aci][a-][g-quic1275]:0.111111,Sacapulteco[i-quv][w-][a-][g-saca1238]:0.222222,Sipacapense[i-qum][w-qum][a-3121][g-sipa1247]:0.222222)Core_Quichean[i-][w-][a-][g-core1251]:0.111111,(Poqomam[i-poc][w-pcm][a-2317-2319][g-poqo1253]:0.222222,Poqomchi`[i-poh][w-][a-2318][g-poqo1254]:0.222222)Pocom[i-pob][w-pkm][a-][g-poco1241]:0.111111)Poqom-Quichean[i-][w-][a-][g-poqo1252]:0.111111,Uspanteco[i-usp][w-][a-][g-uspa1245]:0.444444)Greater_Quichean[i-][w-][a-][g-grea1276]:0.111111)Quichean-Mamean[i-][w-][a-][g-quic1274]:0.111111)Core_Mayan[i-][w-][a-][g-core1254]:0.111111,((Itza[i-itz][w-itz][a-1660][g-itza1241]:0.555556,Mopan_Maya[i-mop][w-mop][a-2058][g-mopa1243]:0.555556)Mopan-Itza[i-][w-][a-][g-mopa1242]:0.111111,((Lacanjá/[i-][w-][a-][g-laca1244]:0.444444,Najá[i-][w-][a-][g-naja1242]:0.444444)Lacandon[i-lac][w-lac][a-1861][g-laca1243]:0.111111,Yucateco[i-yua][w-yct][a-682][g-yuca1254]:0.555556)Yucatec-Lacandon[i-][w-][a-][g-yuca1253]:0.111111)Yucatecan[i-][w-][a-][g-yuca1252]:0.111111)Yucatecan-Core_Mayan[i-][w-][a-][g-yuca1255]:0.111111)Mayan[i-][w-][a-][g-maya1287]:0.111111);", "Mayan");

as.phylo.familytree <- function(tree)
{
  class(tree) <- "phylo";
  tree;
}

# Get name is generic:
get.name <- function(x) UseMethod("get.name");

# Get the tree name:
get.name.familytree <- function( tree )
{
  return (attr(tree, "tree.name"));
}

# Fix the language and group names by replacing special characters with equivalent ASCII characters:
.fix.names <- function( tree, quotes=NA )
{
  if( !inherits(tree, "familytree") ) return (tree);
  # Fix the names:
  .fix.string <- function(s)
  {
    quotes.removed <- FALSE;
    if( !is.na(quotes) && substr(s,1,nchar(quotes))==quotes && substr(s,nchar(s)-nchar(quotes)+1,nchar(s))==quotes ) 
    {
      s <- substr(s,nchar(quotes)+1,nchar(s)-nchar(quotes)); 
      quotes.removed <- TRUE;
    }
    s <- normalize.language.name(s); 
    if( quotes.removed )
    {
      s <- paste0(quotes, s, quotes);
    }
    return (s);
  }
  if( !is.null(attr(tree, "tree.name")) ) attr(tree, "tree.name") <- .fix.string(get.name(tree));
  if( !is.null(tree$tip.label) ) 
    tree$tip.label <- vapply(tree$tip.label, .fix.string, character(1));
  if( !is.null(tree$node.label) ) 
    tree$node.label <- vapply(tree$node.label, .fix.string, character(1));
  return (tree);
}
.fix.names <- cmpfun(.fix.names);

# Convert from the glottolog to my convention using the mapping to find the extra codes:
.convert.glottolog.convention <- function(tree, mapping, quotes="'")
{
  if( !inherits(tree, "familytree") ) return (tree);
  # Convert the tip and node labels (don't touch the root name):
  if( !is.null(tree$tip.label) ) 
    tree$tip.label <- vapply(tree$tip.label, function(s){ .glottolog2mine(s, mapping, quotes) }, character(1));
  if( !is.null(tree$node.label) ) 
    tree$node.label <- vapply(tree$node.label, function(s){ .glottolog2mine(s, mapping, quotes) }, character(1));
  return (tree);
}
.convert.glottolog.convention <- cmpfun(.convert.glottolog.convention);

# Count the number of languages in a family tree:
count.languages <- function( tree )
{
  if( inherits(tree, "familytree") )
  {
    return (Ntip(tree));
  }
  return (0);
}
count.languages <- cmpfun(count.languages);

# Find a rooted tree's root node:
find.root <- function(tree)
{
  if( inherits(tree, "phylo") )
  {
    if( is.rooted(tree) )
    {
      return (Ntip(tree)+1);
    } else
    {
      # Try to apply a quick-and-dirty heuristic to identify the root:
      root.node <- setdiff( 1:(Ntip(tree) + Nnode(tree)), tree$edge[,2] ); # the root has no ancestor
      if( length(root.node) != 1 )
      {
        return (NA);
      } 
      return (root.node);
    }
  }
  return (NA);
}
find.root <- cmpfun(find.root);

# Find a node's parent in a tree or 0 for the root or -1 if something's wrong:
.find.parent <- function(tree, node)
{
  if( inherits(tree, "phylo") && 1 <= node && node <= Nnode(tree, FALSE) )
  {
    parent <- tree$edge[tree$edge[,2] == node,1];
    if( length(parent) == 1 )
    {
      return (parent) 
    } else if( length(parent)==0 )
    {
      return (0);
    } else
    {
      return (-1);
    }
  } else
  {
    return (-1);
  }
}
.find.parent <- cmpfun(.find.parent);

# Get a node's name (or number, if there's no name):
.get.node.name <- function(tree, node)
{
  if( node < 0 || node > Nnode(tree, FALSE) )
  {
    stop("Error retreiving node name!\n");
  } else if(node <= Ntip(tree))
  {
    return (tree$tip.label[node]);
  } else
  {
    if( !is.null(tree$node.label) ) return (tree$node.label[node-Ntip(tree)]) else return (paste0("Node_",node));
  }
}
.get.node.name <- cmpfun(.get.node.name);

# Print a family tree:
print.familytree <- function( tree, as.Newick=FALSE, digits=3, use.ASCII=FALSE, print.brlength=TRUE, show.brlength=TRUE, max.brlength.to.show=100)
{
  if( as.Newick || is.na(root <- find.root(tree)) )
  {
    # Print it as a Newick tree:
    .internal.write.tree <- function(phy, digits = 10) # adapted from ape::.write.tree2() to make sure node names are not altered during writting
    {
      brl <- !is.null(phy$edge.length)
      nodelab <- !is.null(phy$node.label)
      #phy$tip.label <- checkLabel(phy$tip.label)
      #if (nodelab) 
      #  phy$node.label <- checkLabel(phy$node.label)
      f.d <- paste("%.", digits, "g", sep = "")
      cp <- function(x) {
        STRING[k] <<- x
        k <<- k + 1
      }
      add.internal <- function(i) {
        cp("(")
        desc <- kids[[i]]
        for (j in desc) {
          if (j > n) 
            add.internal(j)
          else add.terminal(ind[j])
          if (j != desc[length(desc)]) 
            cp(",")
        }
        cp(")")
        if (nodelab && i > n) 
          cp(phy$node.label[i - n])
        if (brl) {
          cp(":")
          cp(sprintf(f.d, phy$edge.length[ind[i]]))
        }
      }
      add.terminal <- function(i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl) {
          cp(":")
          cp(sprintf(f.d, phy$edge.length[i]))
        }
      }
      n <- length(phy$tip.label)
      parent <- phy$edge[, 1]
      children <- phy$edge[, 2]
      kids <- vector("list", n + phy$Nnode)
      for (i in 1:length(parent)) kids[[parent[i]]] <- c(kids[[parent[i]]], 
                                                         children[i])
      ind <- match(1:max(phy$edge), phy$edge[, 2])
      LS <- 4 * n + 5
      if (brl) 
        LS <- LS + 4 * n
      if (nodelab) 
        LS <- LS + n
      STRING <- character(LS)
      k <- 1
      cp("")
      cp("(")
      getRoot <- function(phy) phy$edge[, 1][!match(phy$edge[, 
                                                             1], phy$edge[, 2], 0)][1]
      root <- getRoot(phy)
      desc <- kids[[root]]
      for (j in desc) {
        if (j > n) 
          add.internal(j)
        else add.terminal(ind[j])
        if (j != desc[length(desc)]) 
          cp(",")
      }
      if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab) 
          cp(phy$node.label[1])
        cp(";")
      }
      else {
        cp(")")
        if (nodelab) 
          cp(phy$node.label[1])
        cp(":")
        cp(sprintf(f.d, phy$root.edge))
        cp(";")
      }
      paste(STRING, collapse = "")
    }
    cat(.internal.write.tree(tree, digits=digits));
  } else
  {
    # Pretty print with indentation:
    if( show.brlength && is.null(tree$edge.length) ) show.brlength <- FALSE; # don't attepmt to print branch length if none is given!
    
    # Auxiliary function to print the indent:
    .show.brlength <- function(brlength.steps, show.spaces=FALSE)
    {
      if( show.brlength && is.numeric(brlength.steps) && brlength.steps > 0 )
      {
        if( brlength.steps < max.brlength.to.show )
        {
          paste0(rep(ifelse(show.spaces," ",ifelse(use.ASCII,"-","\u2500")), brlength.steps), collapse="");
        } else
        {
          if( show.spaces )
          {
            paste0(rep(" ", brlength.steps), collapse="");
          } else
          {
            brlength.steps.left <- floor((brlength.steps - nchar(as.character(brlength.steps))) / 2);
            paste0(paste0(rep(ifelse(use.ASCII,"-","\u2500"), brlength.steps.left), collapse=""), 
                   "/", as.character(brlength.steps), "/", 
                   paste0(rep(ifelse(use.ASCII,"-","\u2500"), brlength.steps - brlength.steps.left - nchar(as.character(brlength.steps))), collapse=""));
          }
        }
      } else
      {
        ifelse(show.spaces," ",ifelse(use.ASCII,"-","\u2500"));
      }
    }
    
    # Recursively print a family tree:
    .print.familytree.recursive <- function(tree, node, digits, indents=NULL, is.last.child=FALSE, use.ASCII=FALSE, print.brlength=TRUE, show.brlength=TRUE, max.brlength.to.show=10)
    {      
      # The branch length printing steps:
      brlength <- NA;
      brlength.steps <- 0;
      
      if( node != find.root(tree) ) # Don't print the root itself!
      {        
        # The branch length:
        if( !is.null(tree$edge.length) && print.brlength )
        {
          the.edge <- (tree$edge[,2] == node); # the edge leading to this node
          if( sum(the.edge) == 0 )
          {
            # The root:
            if( !is.null(tree$root.edge) ) brlength <- tree$root.edge;
          } else if( sum(the.edge) == 1 )
          {
            # Normal node:
            if( !is.null(tree$edge.length) ) brlength <- tree$edge.length[the.edge];
          } else
          {
            stop("Error printing tree!\n");
          }
        }
        
        # Indent:
        brlength.steps <- round(brlength / min.brlength);
        if( !is.null(indents) )
        {
          if( length(indents) > 1 )
          {
            for( i in 1:(length(indents)-1) )
            {
              if( length(grep(ifelse(use.ASCII,"|","\u2502"), indents[i], fixed=TRUE)) > 0 ) cat("  ");
              cat(indents[i]);
            }
          }
          cat(paste0(rep(" ",nchar(indents[length(indents)])),collapse=""));
        }

        if(  !is.last.child )
        {
          cat(paste0(ifelse(use.ASCII,"+","\u251C"),.show.brlength(brlength.steps,FALSE)));
        } else
        {
          cat(paste0(ifelse(use.ASCII,"\\","\u2514"),.show.brlength(brlength.steps,FALSE)));
        }
        
        # The node name:
        cat(.get.node.name(tree, node));
        
        # The branch length:
        if( !is.na(brlength) ) cat(paste0(" : ", sprintf(paste0("%.",digits,"f"),brlength)));
        # End printing the node:
        cat("\n");
      }
      
      # Now go to the children:
      the.children <- which(tree$edge[,1] == node);
      if( length(the.children) > 0 )
      {
        for( i in 1:length(the.children) ) 
          .print.familytree.recursive(tree, 
                                      tree$edge[the.children[i],2], 
                                      digits, 
                                      c(indents, paste0(.show.brlength(brlength.steps-1,TRUE), ifelse(i==length(the.children), "  ", ifelse(use.ASCII,"| ","\u2502 ")))), 
                                      i==length(the.children), 
                                      use.ASCII,
                                      print.brlength,
                                      show.brlength, max.brlength.to.show);
      }
    }
    min.brlength <- NA;
    if( show.brlength )
    {
      # Get the minimum branch length (this will be represented by a single "-"):  
      if( !is.null(tree$edge.length) && sum(tree$edge.length > 0, na.rm=TRUE) > 0 ) min.brlength <- min(tree$edge.length[ tree$edge.length > 0 ], na.rm=TRUE);
    }
    # Print the family name:
    cat(paste0("<",get.name(tree), "> (", count.languages(tree), " tips, ", count.levels(tree)-1, " levels)","\n"));
    # and the its structure:
    .print.familytree.recursive( tree=tree, 
                                 node=find.root(tree), 
                                 digits=digits, 
                                 use.ASCII=use.ASCII, 
                                 print.brlength=print.brlength, 
                                 show.brlength=show.brlength, max.brlength.to.show=max.brlength.to.show );
  }
}
print.familytree <- cmpfun(print.familytree);

as.character.familytree <- function( tree, as.Newick=TRUE, digits=3, ... ) return( capture.output( print.familytree(tree, as.Newick, digits, ...) ) );

# Extract the languages from a family tree:
extract.languages <- function( tree )
{
  if( inherits(tree, "familytree") )
  {
    return (tree$tip.label);
  }
  return (NULL);
}
extract.languages <- cmpfun(extract.languages);

# Extract the internal nodes from a family tree:
extract.internal.nodes <- function( tree )
{
  if( inherits(tree, "familytree") )
  {
    return (tree$node.label);
  }
  return (NULL);
}
extract.internal.nodes <- cmpfun(extract.internal.nodes);

# Count the levels in a family tree:
count.levels <- function( tree )
{
  if( !is.na(find.root(tree)) )
  {
    # Recursivelly count its levels:
    .count.levels <- function( tree, node )
    {
      the.children <- which(tree$edge[,1] == node);
      if( length(the.children) > 0 )
      {
        levels <- vapply(the.children, function(i) .count.levels(tree, tree$edge[i,2]), numeric(1));
        return (max(levels)+1);
      } else
      {
        return (1);
      }
    }
    return (.count.levels(tree, find.root(tree)));
  } else
  {
    return (0);
  }
}
count.levels <- cmpfun(count.levels);

# Extract all subtrees of a given level (1=root) from a family tree and return them as a vector of nodes:
extract.subtrees.of.level <- function( tree, level=2 )
{
  # Sainty check:
  if( count.levels(tree) < level )
  {
    return( NULL );
  } else
  {
    # Recursively extract the subtrees:
    .extract.subtrees.of.level <- function( tree, node, cutoff.level=2, node.level=1 )
    {
      if( cutoff.level == node.level )
      {
        # This is a good subtree:
        x <- NULL; try(x <- extract.clade(tree, node, 1), silent=TRUE);
        return (x);
      } else
      {
        the.children <- which(tree$edge[,1] == node);
        if( length(the.children) > 0 )
        {
          subtrees <- list();
          for( i in the.children )
          {
            x <- .extract.subtrees.of.level(tree, tree$edge[i,2], cutoff.level=cutoff.level, node.level=node.level+1);
            # Flatten the list and remove NULLs:
            if( is.null(x) || length(x) == 0 )
            {
              next;
            } else if( inherits(x, "familytree") )
            {
              subtrees[[length(subtrees)+1]] <- x;
            } else
            {
              for( j in 1:length(x) ) subtrees[[length(subtrees)+1]] <- x[[j]];
            }
          }
          return (subtrees);
        } else
        {
          return (NULL);
        }
      }
    }   
    return (.extract.subtrees.of.level(tree, find.root(tree), cutoff.level=level+1));
  }
}
extract.subtrees.of.level <- cmpfun(extract.subtrees.of.level);

# Prune a tree by keeping only those languages and internal nodes in the given set:
prune.family.to.subset <- function( tree, languages.set )
{ 
  if( inherits(tree, "familytree") && !is.null(languages.set) && length(languages.set) > 0 )
  {
    # Are there internal nodes in this list?
    internal.nodes <- intersect(extract.internal.nodes(tree), languages.set);
    if( length(internal.nodes) > 0 )
    {
      # Extract the paths to all the nodes (internal and terminal):
      nodes.paths <- lapply(languages.set, function(s) extract.path(tree, s, include.root=TRUE, include.brlen=TRUE)); names(nodes.paths) <- languages.set;
    
      # Build a new tree from these paths:
      subtree <- NULL;
      for( path in nodes.paths ) subtree <- add.tree.path(subtree, path[-1], brlens=as.numeric(names(path)[-1]));
    } else
    {
      # Standard tip prunning:
      subtree <- drop.tip(tree, setdiff(tree$tip.label, languages.set));
    }

    return (subtree);
  } else
  {
    return (NULL);
  }
}
prune.family.to.subset <- cmpfun(prune.family.to.subset);

# Extract the path from the root to a node (language or internal node) as a vector of strings starting with the root; if requested, the brlens are included as the names:
extract.path <- function( tree, node, include.root=TRUE, include.brlen=TRUE )
{ 
  if( inherits(tree, "familytree") && length(node)==1 && !is.null(root <- find.root(tree)) )
  {
    # Identify the node:
    if( is.numeric(node) ) {
      cur.node <- node;
    } else if( node %in% tree$tip.label ) { 
      cur.node <- which(tree$tip.label == node);
    } else if( node %in% tree$node.label ) {
      cur.node <- Ntip(tree)+which(tree$node.label == node);
    } else {
      cur.node <- (-1);
    }
    if( cur.node < 0 || cur.node > Nnode(tree, internal.only=FALSE) ) return (NULL);
    # Walk up from the tip to the root:
    path <- NULL; prev.node <- (-1);
    while(TRUE)
    {
      path <- c(.get.node.name(tree, cur.node), path);
      if( include.brlen && !is.null(tree$edge.length) && prev.node != (-1) ) names(path)[2] <- tree$edge.length[ tree$edge[,2] == prev.node & tree$edge[,1] == cur.node ];
      prev.node <- cur.node;
      cur.node <- .find.parent(tree, cur.node);
      if( cur.node <= 0 || (!include.root && cur.node == root) ) break;
    }
    return (path);
  } else
  {
    return (NULL);
  }
}
extract.path <- cmpfun(extract.path);

# Add a path (root -> tip, given as a vector of node names, possibly with branch lengths given as a vector of numbers [the original brlenghts have priority for the shared path]) to an existing tree:
add.tree.path <- function( tree, path, brlens=NULL, warn.on.duplicated.tip=TRUE )
{
  if( is.null(path) || is.na(path) || length(path)==0  )
  {
    # Nothing to add!
    return (tree);
  }
  path.len <- length(path);
  
  if( is.null(tree) )
  {
    # Create a new tree from this path:
    if( path.len == 1 )
    {
      # Just a tip: this is a special case to be treated differently:
      if( !is.null(brlens) )
      {
        new.tip <- list(edge=matrix(c(2,1),1,2),
                        tip.label=path[path.len],
                        edge.length=brlens[path.len],
                        Nnode=1);
      } else
      {
        new.tip <- list(edge=matrix(c(2,1),1,2),
                        tip.label=path[path.len],
                        Nnode=1);
      }
      class(new.tip)<-c("familytree","phylo");
      attr(new.tip, "tree.name") <- path[1];
      return (new.tip);
    } else
    {
      # Add a whole path including at least one internal node:
      newick.path <- rep("(",path.len);
      for( i in path.len:1 )
      {
        newick.path <- c(newick.path, path[i]);
        if( !is.null(brlens) ) newick.path <- c(newick.path, ":", brlens[i]);
        newick.path <- c(newick.path, ")");
      }
      newick.path <- c(newick.path, ";")
      new.branch <- read.newick(text=paste0(newick.path,collapse=""));
      class(new.branch) <- append("familytree", class(new.branch));
      attr(new.branch, "tree.name") <- path[1];
      return (new.branch);
    }
  }
  
  if( !inherits(tree,"familytree") || is.na(find.root(tree)) || is.null(tree$node.label) )
  {
    warning( paste0("Adding path (", paste0(path,collapse=", "), ") to non-tree failed\n") );
    return (NULL);
  }
  
  if( warn.on.duplicated.tip && path[length(path)] %in% tree$tip.label ) warning(paste0("Adding path '", paste0(path,collapse=","), "' to the tree: the tip node '", path[length(path)], "' already present in the tree\n"));
  
  # Find the highest node already in the tree following the path from the root:
  root <- find.root(tree); children <- which(tree$edge[,1] == root); NTotNodes <- Ntip(tree)+Nnode(tree);
  if( !is.null(brlens) && length(brlens)==1 ) brlens <- rep(brlens,path.len);
  if( .get.node.name(tree,root) == "" && length(children) == 1 )
  {
    # The degenrate un-named root with a single descendant: go to this unique descendant:
    cur.node <- tree$edge[children[1],2];
  } else
  {
    # Go to the root itself:
    cur.node <- root;
  }
  cur.path <- 1;
  found <- FALSE;
  while( cur.node >= 1 && cur.node <= NTotNodes && cur.path <= path.len && path[cur.path] == .get.node.name(tree,cur.node) )
  {
    if( cur.node <= Ntip(tree) )
    {
      # This is a tip:
      if( cur.path != path.len ) warning(paste0("Path '", paste0(path,collapse=","), "' contains a tip node (", .get.node.name(tree,cur.node), ")\n" ));
      return (tree); # nothing to add here!
    } else
    {
      # Still an internal node: look for the child that matches the next path element:
      if( cur.path == path.len )
      {
        # No next path element: nothing to add:
        return (tree);
      }
      cur.path <- cur.path+1;
      children <- which(tree$edge[,1] == cur.node);
      if( length(children) == 0 ){ warning("Malformed tree: internal node without children!\n"); return (NULL); }
      matching.children <- vapply(children, function(i) .get.node.name(tree,tree$edge[i,2]) == path[cur.path], logical(1));
      if( sum(matching.children,na.rm=TRUE) >= 1 )
      {
        cur.node <- tree$edge[children[which(matching.children)[1]],2]; # pick the first matching child and continue this search
        next;
      } else
      { 
        # Good, no child is matching this path element: add it to the tree then!
        found <- TRUE;
        break;
      }
    }
  }
  
  if( found )
  {
    # Ok, so we have to add the path starting with cur.path to the node cur.node:
    if( cur.path == path.len )
    {
      # Add just a tip: this is a special case to be treated differently:
      if( !is.null(brlens) )
      {
        new.tip <- list(edge=matrix(c(2,1),1,2),
                        tip.label=path[path.len],
                        edge.length=brlens[path.len],
                        Nnode=1);
      } else
      {
        new.tip <- list(edge=matrix(c(2,1),1,2),
                        tip.label=path[path.len],
                        Nnode=1);
      }
      class(new.tip)<-"phylo";
      new.tree <- bind.tree(tree, new.tip, where=cur.node);
    } else
    {
      # Add a whole path including at least one internal node:
      newick.path <- rep("(",path.len-cur.path+1);
      for( i in path.len:cur.path )
      {
        newick.path <- c(newick.path, path[i]);
        if( !is.null(brlens) ) newick.path <- c(newick.path, ":", brlens[i]);
        newick.path <- c(newick.path, ")");
      }
      newick.path <- c(newick.path, ";")
      new.branch <- read.newick(text=paste0(newick.path,collapse=""));
      new.tree <- bind.tree(tree, new.branch, where=cur.node);
    }
    return (new.tree);
  }
}
add.tree.path <- cmpfun(add.tree.path);
# add.tree.path( tree=familytree(text="(((a:1,b:2)c:1,(d:4,e:5)f:2)g:1);",tree.name="test"), path=c("g","f","h","i"), brlens=c(7,2,10,3) )
# add.tree.path( tree=NULL, path=c("g","f","h","i"), brlens=c(7,2,10,3) )

# Test if two familytree objects are equal:
`==.familytree` <- function(tree1, tree2) all.equal(tree1, tree2, use.edge.lengt=TRUE, use.tip.label=TRUE);

# Collapse single nodes and the reverse (restore the collapsed single nodes)
# Insert a row in a matrix at a given position:
.insert.row <- function(m, row=rep(NA,ncol(m)), where=nrow(m))
{
  if( is.null(m) || !is.matrix(m) || where <= 0 || length(row) != ncol(m) ) return (m); # nothing to do
  rnames <- rownames(m);
  if( where == 1 )
  {
    m <- rbind(row,m);
  } else if( where <= nrow(m) )
  {
    m <- rbind(m[1:(where-1),], row, m[where:nrow(m),]);
  } else if( where == nrow(m)+1 )
  {
    m <- rbind(m,row);
  } else
  {
    m <- rbind(m,matrix(NA,ncol=ncol(m),nrow=where-nrow(m)-1),row);
  }
  if( is.null(rnames) ) rownames(m) <- NULL;
  return (m);
}
.insert.row <- cmpfun(.insert.row);

# Insert an element in a vector at a given position:
.insert.element <- function(v, element=NA, where=length(v))
{
  if( is.null(v) || where <= 0 ) return (v); # nothing to do
  rnames <- names(v);
  if( where == 1 )
  {
    v <- c(element,v);
  } else if( where <= length(v) )
  {
    v <- c(v[1:(where-1)], element, v[where:length(v)]);
  } else if( where == length(v)+1 )
  {
    v <- c(v,element);
  } else
  {
    v <- c(v,rep(NA,where-length(v)-1),element);
  }
  if( is.null(rnames) ) names(v) <- NULL;
  return (v);
}
.insert.element <- cmpfun(.insert.element);

# Collapse single nodes also returning the collapsed nodes so that the original topology can be reconstructed later on:
collapse.singles.reversible <- function(tree) # return a list containing the new tree and the collapsed singles
{
  # This is shamelessly inspired from ape's original collapse.fingles() function
  # Cache several tree properties:
  elen <- tree$edge.length;
  xmat <- tree$edge;
  node.lab <- tree$node.label;
  nnode <- tree$Nnode;
  ntip <- length(tree$tip.label);
  
  # The original tree's root:
  root <- find.root(tree);
  
  # Start processing the singles:
  singles <- NA;
  removed.singles <- data.frame("node"=rep(NA,nnode), "name"=NA, "prev.node"=NA, "next.node"=NA, "prev.brlen"=NA, "next.brlen"=NA); k <- 1; tree.orig <- tree; 
  while( length(singles) > 0 )
  {
    tx <- tabulate( xmat[,1] );
    singles <- which( tx == 1 ); # singles are those nodes with just one descendant
    if( length(singles) > 0 )
    {
      i <- singles[1]; # focus on the first single
      
      prev.node <- which( xmat[,2] == i ); # the ancestor
      next.node <- which( xmat[,1] == i ); # the single descendant
      prev.single.edge <- which(xmat[,1] == i); # the (ancestor -> single) edge

      xmat[ xmat > i ] <- xmat[ xmat > i ] - 1L; # adjust the node numbering

      # Store this to be removed node:
      removed.singles$node[k] <- i;
      if( !is.null(node.lab) ){ removed.singles$name[k] <- node.lab[i - ntip]; };
      if( length(prev.node) > 0 ){ removed.singles$prev.node[k] <- xmat[prev.node,1]; removed.singles$prev.brlen[k] <- elen[prev.node]; }
      if( length(next.node) > 0 ){ removed.singles$next.node[k] <- xmat[next.node,2]; removed.singles$next.brlen[k] <- elen[next.node]; }
      
      # Connect the (ancestor) directly to the (descendant) removing the (single):
      xmat[prev.node, 2] <- xmat[next.node, 2]; # connect the ancestor directly to the descendant        
      xmat <- xmat[ -prev.single.edge,]; # delete the (ancestor -> single) edge

      elen[prev.node] <- elen[prev.node] + elen[next.node]; # the new direct branch's length is the sum of the (ancestor -> single) and the (single -> descendant)'s branch lengths 
      if( !is.null(node.lab) ) node.lab <- node.lab[-c(i - ntip)]; # adjst the node labels as well
      nnode <- nnode - 1L; # one less node
      elen <- elen[-next.node]; # and one less edge
      
      k <- k+1; # prepare for the (possible) next removed single
    }
  }
  removed.singles <- removed.singles[ !is.na(removed.singles$node), ]; # clean them
  
  # Update the tree...
  tree$edge <- xmat;
  tree$edge.length <- elen;
  tree$node.label <- node.lab;
  tree$Nnode <- nnode;
  
  # ...and return it
  list("original.tree"=tree.orig, "collapsed.tree"=tree, "removed.singles"=removed.singles, "original.root"=root);
}
collapse.singles.reversible <- cmpfun(collapse.singles.reversible);

# Reverse the collapsing of single nodes:
reverse.collapse.singles <- function(collapsed.tree, reverse.info, 
                                     restore.brlen.method=c("original.proportion","equal.proportion","first.zero")) # use the original.tree and collapsed.singles list to restore the collapsed.tree
{
  if( is.null(reverse.info) || is.null(reverse.info$removed.singles) || nrow(reverse.info$removed.singles) == 0 ) return (collapsed.tree); # nothing to restore!
  
  restore.brlen.method <- restore.brlen.method[1];
  
  # Cache several tree properties:
  elen <- collapsed.tree$edge.length;
  xmat <- collapsed.tree$edge;
  node.lab <- collapsed.tree$node.label;
  nnode <- collapsed.tree$Nnode;
  ntip <- length(collapsed.tree$tip.label);
  singles <- reverse.info$removed.singles; # less typing
  
  # The removed.singles list records the order of removal: start from the end
  for( i in nrow(singles):1 )
  {
    s <- singles[i,]; # save typing...
    # Try to locate the ancestor ("prev") and decendant ("next") nodes in the tree as we need to restore this single between them:
    prev.node <- s$prev.node; next.node <- s$next.node;
    if( is.na(prev.node) )
    {
      # This is the new root: insert it:
      xmat <- .insert.row(xmat, c(s$node, next.node), 1); # insert the link
      
      # Adjust the node numbering and names:
      xmat[ xmat >= s$node ] <- xmat[ xmat >= s$node ] + 1L; xmat[1,1] <- s$node; # but make sure to keep the right root node
      if( !is.null(node.lab) ) node.lab <- .insert.element(node.lab, s$name, s$node - ntip);
      nnode <- nnode + 1L;
      
      # Update the brlens:
      if( !is.null(elen) )
      {
        elen <- .insert.element(elen, (s$next.brlen), 1);
      }
    } else
    {
      # Regular internal node:
      dlink <- which(xmat[,1] == prev.node & xmat[,2] == next.node); # locate the (ancestor -> descendant) direct link
      if( length(dlink) != 1 )
      {
        warning(paste0("Error restoring collapsed single node ",s$orig.node,ifelse(!is.na(s$node.name),paste0(", '",s$node.name,"'"),"")," cannot locate branch!\n"));
        return (collapsed.tree);
      }
      
      # Adjust the node numbering and names:
      xmat[ xmat >= s$node ] <- xmat[ xmat >= s$node ] + 1L;
      if( !is.null(node.lab) ) node.lab <- .insert.element(node.lab, s$name, s$node - ntip);
      nnode <- nnode + 1L;
      
      # Replace the direct link by two links (ancestor -> single) and (single -> ancestor):
      xmat[dlink,2] <- s$node; 
      # Insert a new link (single -> descendant):
      xmat <- .insert.row(xmat, c(s$node, ifelse(next.node >= s$node, next.node+1L, next.node)), dlink+1);
      
      # Update the brlens:
      if( !is.null(elen) )
      {
        if( !(restore.brlen.method %in% c("original.proportion","equal.proportion","first.zero")) )
        {
          warning("Unknown branch length restoration method: defaulting to original proportion\n");
          restore.brlen.method <- "original.proportion";
        }
        k <- elen[dlink]; # the brlen to be split
        if( restore.brlen.method == "original.proportion" )
        {
          elen[dlink] <- k * s$next.brlen / (s$prev.brlen + s$next.brlen);
        } else if( restore.brlen.method == "equal.proportion" )
        {
          elen[dlink] <- k / 2;
        } else if( restore.brlen.method == "first.zero" )
        {
          elen[dlink] <- 0;
        } else
        {
          # This should never happen
        }
        elen <- .insert.element(elen, (k - elen[dlink]), dlink);
      }
    }
  }
  
  # Update the tree...
  collapsed.tree$edge <- xmat;
  collapsed.tree$edge.length <- elen;
  collapsed.tree$node.label <- node.lab;
  collapsed.tree$Nnode <- nnode;
  
  # ... and return it:
  return (collapsed.tree);
}
reverse.collapse.singles <- cmpfun(reverse.collapse.singles);
# TEST: tree <- familytree("(((a:0.1)A:0.5,(b1:0.2,b2:0.1)B:0.2,(c:0.2)C:0.1)R:0.1);", "test"); # the test tree
# TEST: tree <- familytree("((a:0.1)A:0.5,(b1:0.2,b2:0.1)B:0.2);", "test"); # the test tree
# TEST: tree <- familytree("((((('Maya. Yucatec [i-yua][w-yct][a-682][g-yuca1254]':1,'Lacandon [i-lac][w-lac][a-1861][g-laca1243]':1)'Yucatec-Lacandon [i-][w-][a-][g-]':1,('Maya. Mopan [i-mop][w-mop][a-2058][g-mopa1243]':1,'Itza` [i-itz][w-itz][a-1660][g-itza1241]':1)'Mopan-Itza [i-][w-][a-][g-]':1)'Yucatecan [i-][w-][a-][g-]':1,((('Q`anjob`al [i-kjb][w-kea][a-1760][g-qanj1241]':1,'Jakalteko [i-jac][w-jak][a-460][g-popt1235]':1,'Akateko [i-knj][w-kwe][a-1787][g-west2635]':1)'Q`anjob`al-Akateko-Jakalteko [i-][w-][a-][g-]':1,'Mocho [i-mhc][w-][a-][g-moch1257]':1)'Q`anjob`alan [i-][w-][a-][g-]':1,('Tojolabal [i-toj][w-toj][a-2598][g-tojo1241]':1,'Chuj [i-cac][w--chj][a-1107-][g-chuj1250-chuj1249]':1)'Chujean [i-][w-][a-][g-]':1)'Q`anjob`alan-Chujean [i-][w-][a-][g-]':1,((('Tektiteko [i-ttc][w-tec][a-2627][g-tekt1235]':1,'Mam [i-mam][w-][a-2039-2101][g-mamm1241]':1)'Teco-Mam [i-][w-][a-][g-]':1,('Ixil [i-ixl][w-][a-1665][g-ixil1251]':1,'Awakateko [i-agu][w-agu][a-855][g-agua1252]':1)'Awakateko-Ixil [i-][w-][a-][g-]':1)'Mamean [i-][w-][a-][g-]':1,(((('Poqomchi` [i-poh][w-][a-2318][g-poqo1254]':1)'Poqomchi` [i-][w-][a-][g-]':1,('Poqomam [i-poc][w-pcm][a-2317-2319][g-poqo1253]':1)'Poqomam [i-][w-][a-][g-]':1)'Poqom [i-][w-][a-][g-]':1,(('Tz`utujil [i-tzj][w--tzu][a-388-][g-tzut1248-cakc1244]':1,'Kaqchikel [i-cak][w-][a-][g-kaqc1270]':1)'Kaqchikel-Tz`utujil [i-][w-][a-][g-]':1,'Sipakapense [i-qum][w-qum][a-3121][g-sipa1247]':1,'Sakapulteko [i-quv][w-][a-][g-saca1238]':1,'K`iche` [i-quc][w-qch][a-337][g-kich1262]':1,'Achi [i-acr][w-][a-826][g-achi1256]':1)'Core K`ichean [i-][w-][a-][g-]':1)'Poqom-K`ichean [i-][w-][a-][g-]':1,'Uspanteko [i-usp][w-][a-][g-uspa1245]':1,'Q`eqchi` [i-kek][w-kek][a-1730][g-kekc1242]':1)'K`ichean [i-][w-][a-][g-]':1)'K`ichean-Mamean [i-][w-][a-][g-]':1,(((('Tzotzil [i-tzo][w-][a-2652-2654-2655][g-tzot1259]':1)'Tzotzil [i-][w-][a-][g-]':1,('Tzeltal [i-tzh][w-tza-tzt][a-1433-1434-2651-804][g-tzel1254]':1)'Tzeltal [i-][w-][a-][g-]':1)'Tzeltalan [i-][w-][a-][g-]':1,(('Ch`orti` [i-caa][w-coi][a-1105][g-chor1273]':1)'Chorti-Cholti [i-][w-][a-][g-]':1,(('Chol [i-ctu][w-][a-1196][g-chol1282]':1)'Chol [i-][w-][a-][g-]':1,'Chontal. Tabasco [i-chf][w-cmy][a-1136][g-taba1266]':1)'Chol-Chontal [i-][w-][a-][g-]':1)'Cholan [i-][w-][a-][g-]':1)'Cholan-Tzeltalan [i-][w-][a-][g-]':1)'Core Mayan [i-][w-][a-][g-]':1)'Yucatecan-Core Mayan [i-][w-][a-][g-]':1,(('Huastec [i-hus][w-][a-619][g-huas1242]':1)'Huastec [i-][w-][a-][g-]':1,'Chicomuceltec [i-cob][w-cec][a-1167][g-chic1271]':1)'Huastecan [i-][w-][a-][g-]':1)'Mayan [i-][w-][a-][g-]':1);", "test"); # a really complex test tree


# Compare (compute the distance) between two trees using the following distances:
#   - PH85 as implemented by ape's dist.topo;
#   - score as implemented by ape's dist.topo
compare.trees <- function( tree1, tree2 )
{
  # The return value:
  d <- c("PH85"=NA,"score"=NA);
  
  # Preconditions:
  if( is.null(tree1) || is.null(tree2) || !inherits(tree1,"familytree") || !inherits(tree2,"familytree") )
  {
    return (d);
  }
  
  # PH85 and score:
  try( d["PH85"]  <- abs(dist.topo( tree1, tree2, method="PH85"  )), silent=TRUE );
  try( d["score"] <- abs(dist.topo( tree1, tree2, method="score" )), silent=TRUE );
  
  return (d);
}
compare.trees <- cmpfun(compare.trees);

# Compute and apply various types of branch length methods.
# Return the tree (possibly NULL) and comments/supplementary info.
# Methods:
#   - none: all branch length are NA
#   - constant: all root-to-tip lengths are equal
#   - proportional: all branches have the same leagth so that the root-to-tip path is proportional to number of splits
#   - grafen: grafen's method
#   - nnls: apply a given tip-to-tip distances matrixes with minimal distortion; if it has rown- & colnames these must be the fully described language names 'name [i-][w-][g-]', otherwise they are assumed to match the languages as given by the preorder exploration of the tree
#   - ga: as nnls but using a genetic algorithm approach to minimizing the error
#   - nj: build the neighbour-joining tree of the languages with the given distances matrix (the internal nodes do not have names in this case)
# GA can set a branch length to NA (when not enough info is in the distances matrix to actually estimate it) but some programs don't like NA branches; if so replace NA length by replace.NA.brances.with
# Standard phylogenetic methods such as nnls must collapse single nodes, but setting restore.collapsed.singles to TRUE makes sure they are restored afterwards
compute.brlen <- function( tree, method=c("none","constant","proportional","grafen","nnls","ga","nj"), constant=1.0, distmatrix=NULL, replace.NA.brlen.with=NA, restore.collapsed.singles=TRUE )
{
  # Checks:
  if( !inherits(tree, "familytree") )
  {
    return (list("root"=NULL,"comment"="NULL or non-familytree object"));
  }
  
  # Auxiliary specialized methods:
  
  # No branch lengths:
  .compute.brlen.none <- function( tree )
  {
    tree$edge.length <- NULL;
    return (list("tree"=tree,"comment"="No branch lengths"));
  }

  # Total path legth proportional to number of splits (== all branch have same length):
  .compute.brlen.proportional <- function( tree, constant, replace.NA.brlen.with=NA )
  {
    new.tree <- ape::compute.brlen(tree, constant);
    
    # Replace any NA branch length by the given value (if any):
    if( !is.na(replace.NA.brlen.with) ){ new.tree$edge.length[ is.na(new.tree$edge.length) ] <- replace.NA.brlen.with; } # replace NA branch length by the value requested

    return (list("tree"=new.tree,"comment"=paste("Path length proportional to the number of splits with atomic branch length = ",constant,sep="")));
  }
  
  # Total path legth is constant (== the same amount of evolution on all branches):
  .compute.brlen.constant <- function( tree, constant, replace.NA.brlen.with=NA )
  {
    # Checks:
    if( constant <= 0.0 )
    {
      return (list("tree"=NULL,"comment"=paste("Path length ",constant," must be a positive number!",sep="")));
    }
    
    # First compute the number of branches on each path:
    tmp <- .compute.brlen.proportional( tree, 1.0 );
    if( is.null(tmp$tree) )
    {
      return (tmp);
    }
    tree <- tmp$tree;
    
    # Get the tree depth:
    tree.levels <- count.levels( tree );
    
    # Get the smallest step:
    brlen <- constant / tree.levels;
    
    .compute.brlen.constant.recursive <- function( node, brlen, remaning.path.len )
    {
      if( remaning.path.len > 0.0 && length(the.children <- which(tree$edge[,1] == node)) == 0 )
      {
        # Language node:
        tree$edge.length[tree$edge[,2] == node] <<- remaning.path.len;
      } else
      {
        # Visit the children:
        tree$edge.length[tree$edge[,2] == node] <<- brlen;
        for( i in the.children )
        {
          .compute.brlen.constant.recursive(tree$edge[i,2],  brlen, remaning.path.len - brlen);
        }
      }
    }
    
    .compute.brlen.constant.recursive(find.root(tree), brlen, constant+brlen );
        
    # Replace any NA branch length by the given value (if any):
    if( !is.na(replace.NA.brlen.with) ){ tree$edge.length[ is.na(tree$edge.length) ] <- replace.NA.brlen.with; } # replace NA branch length by the value requested

    return (list("tree"=tree, "comment"=paste("Path length is constant ",constant,sep="")));
  }
  
  # Grafen's (1989) method: each node is given a ‘height’, namely the number of leaves of the subtree minus one, 0 for leaves. 
  # Branch lengths are then computed as the difference between height of lower node and height of upper node.
  .compute.brlen.grafen <- function( tree, replace.NA.brlen.with=NA )
  {
    # Unfortunately ape::compute.brlen(tree, "Grafen") apparently cannot deal with certain tree topologies (e.g., root -> node1 -> (leaf1, leaf2)) :(
    # return (list("tree"=ape::compute.brlen(tree, "Grafen"), "comment"="Grafen's (1989) method"));
    
    # So I am reimplementing it here:
    tree <- .compute.brlen.proportional(tree, 1.0)$tree; # make sure we have the branch length structure initialized
    
    # First, compute the heights:
    heights <- rep(NA, Nnode(tree, FALSE));
    .grafen.compute.node.heights <- function( tree, node )
    {
      if( node <= Ntip(tree) )
      {
        heights[node] <<- 0;
        return (0);
      } else
      {
        height <- 0;
        the.children <- which(tree$edge[,1] == node);
        if( length(the.children) > 0 )
        {
          for( i in 1:length(the.children) ) 
          {
            height <- height + .grafen.compute.node.heights( tree, tree$edge[the.children[i],2] );
          }
        }
        heights[node] <<- height + length(the.children) - 1;
        return (heights[node]);
      }
    }
    .grafen.compute.node.heights(tree, find.root(tree) );
    
    # Then compute the branch lengt by substracting the height of the lower and upper nodes:
    for( i in 1:nrow(tree$edge) )
    {
      tree$edge.length[i] <- heights[ tree$edge[i,1] ] - heights[ tree$edge[i,2] ];
    }
        
    # Replace any NA branch length by the given value (if any):
    if( !is.na(replace.NA.brlen.with) ){ tree$edge.length[ is.na(tree$edge.length) ] <- replace.NA.brlen.with; } # replace NA branch length by the value requested
    
    return (list("tree"=tree,"comment"="Grafen's (1989) method"));
  }
  
  # Map a given distances matrix between languages to the tree with minimal distortion (implementation based on nnls.tree in package phangorn):
  .compute.brlen.nnls <- function( tree, distmatrix, replace.NA.brlen.with=NA, restore.collapsed.singles=TRUE )
  {
    # Check the matrix dimensions and match to the languages:
    lgs <- extract.languages( tree );
    if( is.null(lgs) || length(lgs) <= 1 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "Too few languages...\n" );
      return (list("tree"=NULL,"comment"=paste("NNLS: too few languages",sep="")));
    }
    
    if( nrow(distmatrix) != ncol(distmatrix) )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix must be square!\n" );
      return (list("tree"=NULL,"comment"=paste("NNLS: the distances matrix must be square",sep="")));
    }
    if( is.null(rownames(distmatrix)) && is.null(colnames(distmatrix)) )
    {
      if( nrow(distmatrix) != length(lgs) )
      {
        if( PRINT_DEBUG_MESSAGES ) cat( "A distances matrix without rown and column names must has the same number of languages as those in the tree!\n" );
        return (list("tree"=NULL,"comment"=paste("NNLS: a distances matrix without rown and column names must have the same number of languages as those in the tree",sep="")));
      } else
      {
        # Assume the language order in the matrix is the preorder in the tree:
        rownames(distmatrix) <- colnames(distmatrix) <- lgs;
      }
    } else
    {
      # Only one is absent: fill it rom the other one:
      if( is.null(rownames(distmatrix)) )
      {
        rownames(distmatrix) <- colnames(distmatrix);
      } else if( is.null(colnames(distmatrix)) )
      {
        colnames(distmatrix) <- rownames(distmatrix);
      }
    }
    
    # Extract the submatrix corresponding to these languages:
    distmatrix <- distmatrix[ rownames(distmatrix) %in% lgs, colnames(distmatrix) %in% lgs ];
    if( is.null(distmatrix) || sum(!is.na(distmatrix)) == 0 || class(distmatrix) != "matrix"  || nrow(distmatrix) <= 1 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix is too small...\n" );
      return (list("tree"=NULL,"comment"=paste("NNLS: the distances matrix is too small",sep="")));
    }
    # And order them alphabetically:
    distmatrix <- distmatrix[ order(rownames(distmatrix)), order(colnames(distmatrix)) ];
    # Check that the rows and columns do indeed refer to the same languages:    
    if( sum( rownames(distmatrix) != colnames(distmatrix), na.rm=TRUE ) > 0 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix must refer to the same languages on the rows and columns!\n" );
      return (list("tree"=NULL,"comment"=paste("NNLS: the distances matrix must refer to the same languages on the rows and columns",sep="")));
    }
    
    # Extract the subtree of languages for which we actually have distances:
    subtree <- prune.family.to.subset(tree, rownames(distmatrix));
    # Make sure it has branch lenghts:
    subtree <- .compute.brlen.proportional(subtree, 1.0)$tree;
    if( is.null(subtree) || count.languages(subtree) <= 1 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The selected subtree has too few languages...\n" );
      return (list("tree"=NULL,"comment"=paste("NNLS: the selected subtree has too few languages",sep="")));
    }
    
    # Use nnls.tree in phangorn:
    collapsed.singles.info <- collapse.singles.reversible( subtree ); # collapse singles but make sure we can restore them later
    rescaledphylo <- NULL; nnls.message <- capture.output( try( rescaledphylo <- nnls.tree( distmatrix, collapsed.singles.info$collapsed.tree, rooted=TRUE ), silent=TRUE ) );    
    if( is.null(rescaledphylo) )
    {
      # Some failure, hopefully the message is already printed:
      if( PRINT_DEBUG_MESSAGES ) ( "Error calling nnls.tree()...\n" );
      return (list("tree"=NULL,"comment"=paste("NNLS: error in nnls.tree()",sep="")));
    }
    # Adjust the branch length (sometimes there are small negative values):
    minbrlen <- min( rescaledphylo$edge.length, na.rm=TRUE );
    rescaledphylo$edge.length <- rescaledphylo$edge.length + ifelse( minbrlen < 0.0, -minbrlen, 0.0 );
    
    ## Create the familytree object:
    #class(rescaledphylo) <- append("familytree", class(rescaledphylo));
    #attr(rescaledphylo, "tree.name") <- get.name(tree);
        
    # Replace any NA branch length by the given value (if any):
    if( !is.na(replace.NA.brlen.with) ){ rescaledphylo$edge.length[ is.na(rescaledphylo$edge.length) ] <- replace.NA.brlen.with; } # replace NA branch length by the value requested
    
    # Restore the collapsed singles:
    if( restore.collapsed.singles ) rescaledphylo <- reverse.collapse.singles( rescaledphylo, collapsed.singles.info );
    # DEBUG
    if( !all.equal.phylo(rescaledphylo, subtree, use.edge.length=FALSE) ) stop("reverse.collapse.singles() destroyed the topology for ", get.name(tree),"\n");
    # END DEBUG
    
    return (list("tree"=rescaledphylo,"comment"=paste("NNLS: ",nnls.message,sep="")));
  } 
  
  # Construct the NJ tree that fits the distances matrix:
  .compute.brlen.nj <- function( tree, distmatrix, remove.NA=FALSE, replace.NA.brlen.with=NA )
  {
    # Check the matrix dimensions and match to the languages:
    lgs <- extract.languages( tree );
    if( is.null(lgs) || length(lgs) <= 1 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "Too few languages...\n" );
      return (list("tree"=NULL,"comment"="NJ: too few languages"));
    }
    
    if( nrow(distmatrix) != ncol(distmatrix) )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix must be square!\n" );
      return (list("tree"=NULL,"comment"="NJ: the distances matrix must be square"));
    }
    if( is.null(rownames(distmatrix)) && is.null(colnames(distmatrix)) )
    {
      if( nrow(distmatrix) != length(lgs) )
      {
        if( PRINT_DEBUG_MESSAGES ) cat( "A distances matrix without rown and column names must has the same number of languages as those in the tree!\n" );
        return (list("tree"=NULL,"comment"="NJ: a distances matrix without rown and column names must has the same number of languages as those in the tree"));
      } else
      {
        # Assume the language order in the matrix is the preorder in the tree:
        rownames(distmatrix) <- colnames(distmatrix) <- lgs;
      }
    } else
    {
      # Only one is absent: fill it rom the other one:
      if( is.null(rownames(distmatrix)) )
      {
        rownames(distmatrix) <- colnames(distmatrix);
      } else if( is.null(colnames(distmatrix)) )
      {
        colnames(distmatrix) <- rownames(distmatrix);
      }
    }
    
    # Extract the submatrix corresponding to these languages:
    distmatrix <- distmatrix[ rownames(distmatrix) %in% lgs, colnames(distmatrix) %in% lgs ];
    if( is.null(distmatrix) || sum(!is.na(distmatrix)) == 0 || class(distmatrix) != "matrix"  || nrow(distmatrix) <= 2 ) # NJ requires at least 3 tips in the tree
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix is too small...\n" );
      return (list("tree"=NULL,"comment"="NJ: the distances matrix is too small"));
    }
    if( remove.NA )
    {
      # Keep only the non-NA cells (nj is really sensitive to missing data!):
      distmatrix <- distmatrix[ rowSums(is.na(distmatrix)) == 0, colSums(is.na(distmatrix)) == 0 ];
      if( is.null(distmatrix) || sum(!is.na(distmatrix)) == 0 || class(distmatrix) != "matrix"  || nrow(distmatrix) <= 2 ) # NJ requires at least 3 tips in the tree
      {
        if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix is too small...\n" );
        return (list("tree"=NULL,"comment"="NJ: the distances matrix has too few non-missing data"));
      }
    } else
    {
      tmp <- distmatrix[ rowSums(is.na(distmatrix)) == 0, colSums(is.na(distmatrix)) == 0 ];
      if( is.null(tmp) || sum(!is.na(tmp)) == 0 || class(tmp) != "matrix"  || nrow(tmp) <= 2 ) # NJ requires at least 3 tips in the tree
      {
        if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix is too small...\n" );
        return (list("tree"=NULL,"comment"="NJ: the distances matrix has too few non-missing data"));
      }
    }
    # And order them alphabetically:
    distmatrix <- distmatrix[ order(rownames(distmatrix)), order(colnames(distmatrix)) ];
    # Check that the rows and columns do indeed refer to the same languages:    
    if( sum( rownames(distmatrix) != colnames(distmatrix), na.rm=TRUE ) > 0 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix must refer to the same languages on the rows and columns!\n" );
      return (list("tree"=NULL,"comment"="NJ: the distances matrix must refer to the same languages on the rows and columns"));
    }
    
    # Apply the NJ algorithm to build the tree:
    njphylo <- NULL; njphylo <- njs( distmatrix ); #try( njphylo <- upgma( distmatrix ), silent=TRUE );
    if( is.null(njphylo) )
    {
      # Some failure, hopefully the message is already printed:
      if( PRINT_DEBUG_MESSAGES ) cat( "Error calling nj()...\n" );
      return (list("tree"=NULL,"comment"="NJ: error calling nj()"));
    }
    # Adjust the branch length (sometimes there are small negavitve values):
    minbrlen <- min( njphylo$edge.length, na.rm=TRUE );
    njphylo$edge.length <- njphylo$edge.length + ifelse( minbrlen < 0.0, -minbrlen, 0.0 );
    
    # Create the familytree object:
    class(njphylo) <- append("familytree", class(njphylo));
    attr(njphylo, "tree.name") <- get.name(tree);
        
    # Replace any NA branch length by the given value (if any):
    if( !is.na(replace.NA.brlen.with) ){ njphylo$edge.length[ is.na(njphylo$edge.length) ] <- replace.NA.brlen.with; } # replace NA branch length by the value requested
    
    return (list("tree"=njphylo,"comment"="NJ"));
  } 
  
  # Map a given distances matrix between languages to the tree with minimal distortion (implementation based on genetic algorithms):
  .compute.brlen.ga <- function( tree, distmatrix, method="SSE", replace.NA.brlen.with=NA, restore.collapsed.singles=TRUE )
  {
    # Check the matrix dimensions and match to the languages:
    lgs <- extract.languages( tree );
    if( is.null(lgs) || length(lgs) <= 1 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "Too few languages...\n" );
      return (list("tree"=NULL,"comment"=paste("GA: too few languages",sep="")));
    }
    
    if( nrow(distmatrix) != ncol(distmatrix) )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix must be square!\n" );
      return (list("tree"=NULL,"comment"=paste("GA: the distances matrix must be square",sep="")));
    }
    if( is.null(rownames(distmatrix)) && is.null(colnames(distmatrix)) )
    {
      if( nrow(distmatrix) != length(lgs) )
      {
        if( PRINT_DEBUG_MESSAGES ) cat( "A distances matrix without rown and column names must has the same number of languages as those in the tree!\n" );
        return (list("tree"=NULL,"comment"=paste("GA: a distances matrix without rown and column names must has the same number of languages as those in the tree",sep="")));
      } else
      {
        # Assume the language order in the matrix is the preorder in the tree:
        rownames(distmatrix) <- colnames(distmatrix) <- lgs;
      }
    } else
    {
      # Only one is absent: fill it rom the other one:
      if( is.null(rownames(distmatrix)) )
      {
        rownames(distmatrix) <- colnames(distmatrix);
      } else if( is.null(colnames(distmatrix)) )
      {
        colnames(distmatrix) <- rownames(distmatrix);
      }
    }
    
    # Extract the submatrix corresponding to these languages:
    distmatrix <- distmatrix[ rownames(distmatrix) %in% lgs, colnames(distmatrix) %in% lgs ];
    if( is.null(distmatrix) || sum(!is.na(distmatrix)) == 0 || class(distmatrix) != "matrix"  || nrow(distmatrix) <= 1 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix is too small...\n" );
      return (list("tree"=NULL,"comment"=paste("GA: the distances matrix is too small",sep="")));
    }
    # And order them alphabetically:
    distmatrix <- distmatrix[ order(rownames(distmatrix)), order(colnames(distmatrix)) ];
    # Check that the rows and columns do indeed refer to the same languages:    
    if( sum( rownames(distmatrix) != colnames(distmatrix), na.rm=TRUE ) > 0 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The distances matrix must refer to the same languages on the rows and columns!\n" );
      return (list("tree"=NULL,"comment"=paste("GA: the distances matrix must refer to the same languages on the rows and columns",sep="")));
    }
    
    # Extract the subtree of languages for which we actually have distances:
    subtree <- prune.family.to.subset( tree, rownames(distmatrix) );
    # Make sure it has branch lenghts:
    subtree <- .compute.brlen.proportional( subtree, 1.0 )$tree;
    if( is.null(subtree) || count.languages(subtree) <= 1 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The subtree has too few languages...\n" );
      return (list("tree"=NULL,"comment"=paste("GA: the subtree has too few languages",sep="")));
    }
    
    # Use genetic algorithms:
    # To speed things up, use phylo objects:
    subphylo <- subtree; 
    collapsed.singles.info <- collapse.singles.reversible( subphylo ); # make sure we can restore the collapsed singles later
    subphylo <- collapsed.singles.info$collapsed.tree;
    subphylo <- reorder( subphylo ); # reorder it now to optimize for later (required by cophenetic())
    # How many parameters (branch lengths):
    nparams <- length(subphylo$edge.length);
    if( nparams < 1 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "The subtree has too few branches...\n" );
      return (list("tree"=NULL,"comment"=paste("GA: the subtree has too few branches",sep="")));
    }
    
    # Pre-arrange the distances matrix to fit the cophenetic matrix:
    cpmatrix <- cophenetic.phylo( subphylo );
    if( nrow(distmatrix) != nrow(cpmatrix) || 
        ncol(distmatrix) != ncol(cpmatrix) || 
        !all( sort(rownames(distmatrix)) == sort(rownames(cpmatrix)), na.rm=TRUE ) ||
        !all( sort(colnames(distmatrix)) == sort(colnames(cpmatrix)), na.rm=TRUE ) )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "Error computing the cophenetic distances...\n" );
      return (list("tree"=NULL,"comment"=paste("GA: error computing the cophenetic distances",sep="")));
    }
    distmatrix <- distmatrix[ rownames(cpmatrix), colnames(cpmatrix) ];
    
    # The fitness function compares the cophenetic distances to the actual distances; the methods are like in norm{base} + Squared Sum of Errors (SSE) - optimize by having them separate and compiled!:
    # Optimized cophenetic.phylo: assume the ordering was already done:
    .dist.nodes.optimized <- function(x) 
    {
        n <- Ntip(x);
        m <- x$Nnode;
        nm <- n + m;
        d <- .C(dist_nodes, as.integer(n), as.integer(m), as.integer(x$edge[,1] - 1L), as.integer(x$edge[, 2] - 1L), as.double(x$edge.length), as.integer(Nedge(x)), double(nm * nm), NAOK = TRUE)[[7]];
        dim(d) <- c(nm, nm);
        dimnames(d) <- list(1:nm, 1:nm);
        d;
    }
    .dist.nodes.optimized <- cmpfun(.dist.nodes.optimized);
    .cophenetic.phylo.optimized <- function(x) 
    {
        n <- length(x$tip.label);
        ans <- .dist.nodes.optimized(x)[1:n, 1:n];
        dimnames(ans)[1:2] <- list(x$tip.label);
        ans;
    }
    .cophenetic.phylo.optimized <- cmpfun(.cophenetic.phylo.optimized);
    .fitness.function.norm <- function( x, subphylo, distmatrix, method=c("O", "I", "F", "M", "2") )
    {
      # Apply the new lengths to the phylogeny:
      subphylo$edge.length <- as.numeric(x);
      # Compute the cophenetic distance:
      cpmatrix <- .cophenetic.phylo.optimized( subphylo );
      # Compute the difference to the actual distances matrix:
      diffmatrix <- cpmatrix - distmatrix;
      # And derive the norm:
      return ( -norm( diffmatrix, method ) ); # we want to *reduce* the difference!
    }
    .fitness.function.SSE <- function( x, subphylo, distmatrix, method=c("SSE") )
    {
      # Apply the new lengths to the phylogeny:
      subphylo$edge.length <- as.numeric(x);
      # Compute the cophenetic distance:
      cpmatrix <- .cophenetic.phylo.optimized( subphylo );
      # Compute the difference to the actual distances matrix:
      diffmatrix <- as.numeric( cpmatrix - distmatrix );
      if( all(is.na(diffmatrix)) ) return (Inf) else return ( -sum( diffmatrix^2, na.rm=TRUE ) / length(diffmatrix) ); # we want to *reduce* the difference!
    }
    # Make the choice here and compile it:
    if( method == "SSE" )
    {
      .fitness.function <- cmpfun(.fitness.function.SSE);
    } else
    {
      .fitness.function <- cmpfun(.fitness.function.norm);
    }
    
    min.brlengths <- rep(0.0,nparams); max.brlengths <- rep(max(as.numeric(distmatrix),na.rm=TRUE),nparams);
    # Sanity check:
    if( any(is.na(max.brlengths)) || all(min.brlengths >= max.brlengths) )
    {
      # No range to speak of!
      if( PRINT_DEBUG_MESSAGES ) cat( "The minimum and maximum possible branch length are equal...\n" );
      return (list("tree"=NULL,"comment"=paste("GA: no legal range for branch lengths",sep="")));
    }
    
    maxiter <- GA.MAXITER;
    popSize <- GA.POPSIZE;
    run <- GA.CONSTANTRUN;
    ga.trees <- ga( type = "real-valued", 
                    fitness = .fitness.function,
                    min = rep(0.0,nparams), max = rep(max(as.numeric(distmatrix),na.rm=TRUE),nparams), # individual branches can't get below 0.0 or over the maximum distance between tips
                    popSize = popSize, maxiter = maxiter, run=run, monitor=NULL,
                    subphylo=subphylo, distmatrix=distmatrix, method=method );
    if( is.null(ga.trees) )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "Error doing the GA...\n" );
      return (list("tree"=NULL,"comment"=paste("GA: error doing the genetic algorithm",sep="")));
    }
    # Extract the relevant info:
    if( PRINT_DEBUG_MESSAGES ) cat( "After ", ga.trees@iter, " iterations the error is ", abs(ga.trees@fitnessValue), ".\n" );
    # Build the optimal tree:
    rescaledphylo <- subphylo; rescaledphylo$edge.length <- ga.trees@solution;
    if( any(is.na(rescaledphylo$edge.length)) && PRINT_DEBUG_MESSAGES ){ cat("NA in brlen for family ", get.name(tree), "\n"); }
    
    # Create the familytree object:
    class(rescaledphylo) <- append("familytree", class(rescaledphylo));
    attr(rescaledphylo, "tree.name") <- get.name(tree);

    # Replace any NA branch length by the given value (if any):
    if( !is.na(replace.NA.brlen.with) ){ rescaledphylo$edge.length[ is.na(rescaledphylo$edge.length) ] <- replace.NA.brlen.with; } # replace NA branch length by the value requested
    
    if( restore.collapsed.singles ) rescaledphylo <- reverse.collapse.singles( rescaledphylo, collapsed.singles.info );
    
    # Store the number of iterations as well:
    return (list("tree"=rescaledphylo,"comment"=paste("GA: converged after ",ga.trees@iter," iterations out of ",maxiter," with residual ",ifelse(method=="SSE","SSE",paste(method," norm",sep=""))," error ",abs(ga.trees@fitnessValue),sep="")));
  } 

  
  if( method == "none" )
  {
    return (.compute.brlen.none( tree ));
  } else if( method == "constant" )
  {
    if( is.null(constant) || !is.numeric(constant) || constant < 0.0 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "For method='constant' the constant must be a defined positive number...\n" );
      return (list("tree"=NULL,"comment"="For method='constant' the constant must be a defined positive number"));
    } else
    {
      return (.compute.brlen.constant( tree, constant, replace.NA.brlen.with=replace.NA.brlen.with ));
    }
  } else if( method == "proportional" )
  {
    if( is.null(constant) || !is.numeric(constant) || constant < 0.0 )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "For method='proportional' the constant must be a defined positive number...\n" );
      return (list("tree"=NULL,"comment"="For method='proportional' the constant must be a defined positive number"));
    } else
    {
      return (.compute.brlen.proportional( tree, constant, replace.NA.brlen.with=replace.NA.brlen.with ));
    }
  } else if( method == "grafen" )
  {
    return (.compute.brlen.grafen( tree ));
  } else if( method == "nnls" )
  {
    if( is.null(distmatrix) )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "For method='nnls' the distances matrix must be given...\n" );
      return (list("tree"=NULL,"comment"="For method='nnls' the distances matrix must be given"));
    } else
    {
      return (.compute.brlen.nnls( tree, distmatrix, replace.NA.brlen.with=replace.NA.brlen.with, restore.collapsed.singles=restore.collapsed.singles ));
    }
  } else if( method == "ga" )
  {
    if( is.null(distmatrix) )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "For method='ga' the distances matrix must be given...\n" );
      return (list("tree"=NULL,"comment"="For method='ga' the distances matrix must be given"));
    } else
    {
      return (.compute.brlen.ga( tree, distmatrix, replace.NA.brlen.with=replace.NA.brlen.with, restore.collapsed.singles=restore.collapsed.singles ));
    }
  } else if( method == "nj" )
  {
    if( is.null(distmatrix) )
    {
      if( PRINT_DEBUG_MESSAGES ) cat( "For method='nj' the distances matrix must be given...\n" );
      return (list("tree"=NULL,"comment"="For method='nj' the distances matrix must be given"));
    } else
    {
      return (.compute.brlen.nj( tree, distmatrix, replace.NA.brlen.with=replace.NA.brlen.with ));
    }
  } else
  {
    if( PRINT_DEBUG_MESSAGES ) cat( "Unknown method ", method, "\n" );
    return (list("tree"=NULL,"comment"=paste("Unknown branch length computation method ",method,sep="")));
  }
}
compute.brlen <- cmpfun(compute.brlen);


# Create a collection of family trees (a classification) by:
#   - converting a list of familytree objects, or
#   - reading a Nexus file containing a "trees" block (optionally with a "translate" maps, or
#   - reading a CSV file with at least two columns, one containing the family name and the other the family as a Newick tree (or the string "NULL" or missing data if absent) either by name or number
#         (more parameters will be passed to read.table()).
# The classification can have a name.
languageclassification <- function(classif.name="", familytree.list=NULL, nexus.file="", csv.file="", csv.name.column="Family", csv.tree.column="Tree", sep="\t", quote="", header=TRUE, ...)
{
  if( is.null(familytree.list) && nexus.file == "" && csv.file == "" )
  {
    # Build an empty classification:
    classif <- list();
  } else if( !is.null(familytree.list) )
  {
    # Build it from the list:
    classif <- familytree.list;
  } else if( nexus.file != "" )
  {
    # Parse the Nexus file!
    cat("PARSE NEXUS!\n");
    classif <- list();
  } else
  {
    # Parse the CSV file:
    f <- read.table(file=csv.file, header=header, sep=sep, quote=quote, stringsAsFactors=FALSE, ...);
    if( is.null(f) || nrow(f)==0 || ncol(f)==0 )
    {
      warning("Cannot construct 'languageclassification' from an empty CSV file!\n");
      return (NULL);
    }
    # Get the important columns:
    if( is.numeric(csv.name.column) )
    {
      c.name <- csv.name.column;
    } else
    {
      c.name <- grep(csv.name.column, names(f), fixed=TRUE);
    }
    if( is.numeric(csv.tree.column) )
    {
      c.tree <- csv.tree.column;
    } else
    {
      c.tree <- grep(csv.tree.column, names(f), fixed=TRUE);
    }
    if( length(c.name)!=1 || !is.integer(c.name) || c.name < 0 || c.name > ncol(f) )
    {
      warning("Cannot find the 'name' column in the 'languageclassification' CSV file!\n");
      return (NULL);
    }
    if( length(c.tree)!=1 || !is.integer(c.tree) || c.tree < 0 || c.tree > ncol(f) )
    {
      warning("Cannot find the 'tree' column in the 'languageclassification' CSV file!\n");
      return (NULL);
    }
    # Read the columns and build the list:
    classif <- lapply(1:nrow(f), function(i)
      {
        if( !is.na(f[i,c.tree]) && f[i,c.tree] != "" )
        {
          return (familytree(f[i,c.tree], f[i,c.name]));
        } else
        {
          return (NULL);
        }
      } );
    classif <- classif[!vapply(classif,is.null,logical(1))];
  }
  # The actual object:
  attr(classif, "classif.name") <- classif.name;
  class(classif) <- append("languageclassification", class(classif));
  return (classif);
}
# classif <- languageclassification(classif.name="WALS k=1.0", csv.file="../output/wals/wals-newick-constant=1.00.csv", csv.name.column="Family", csv.tree.column="Tree");

# Get the classification name:
get.name.languageclassification <- function( classif )
{
  return (attr(classif, "classif.name"));
}

print.languageclassification <- function(classif)
{
  if( !inherits(classif, "languageclassification") ) print(classif);
  
  if( length(classif) > 0 )
  {
    cat(paste0("Classification '", get.name(classif), "' contains the following ", length(classif), " families:\n\n", paste0(get.family.names(classif),collapse=", "), "\n\nwith the following structure:\n\n"));
    for( i in seq_along(classif) ){ print(classif[[i]]); cat("\n"); }
  } else
  {
    cat(paste0("Classification '", get.name(classif), "' is empty\n"));
  }
}

# Get all familiy names:
get.family.names <- function( classif )
{
  if( !inherits(classif, "languageclassification") )
  {
    return ( NULL );
  } else
  {
    return (vapply(classif, function(p){ ifelse(is.null(p) || !inherits(p,"familytree"), "", get.name(p)) }, character(1)));
  }
}
get.family.names <- cmpfun(get.family.names);

# Locate a language family in a collection of families:
find.language.family <- function( classif, name )
{
  if( !inherits(classif, "languageclassification") )
  {
    return (NULL);
  } else
  {
    for( i in 1:length(classif) )
    {
      if( get.name(classif[[i]]) == name )
      {
        return (classif[[i]]);
      }
    }
    return (NULL);
  }
}
find.language.family <- cmpfun(find.language.family);

# Add a new language (given as a path from family to language) to a languageclassification:
add.language.to.languageclassification <- function( classif, path, brlens=NULL, warn.on.duplicated.tip=TRUE )
{
  if( is.null(path) || is.na(path) || length(path)==0 )
  {
    cat( "Adding path '", paste0(path,collapse=","), "' to collection failed\n" );
    return (NULL);
  }
  path.len <- length(path);
  if( !is.null(brlens) && length(brlens)==1 ) brlens <- rep(brlens, path.len);
  
  if( is.null(classif) )
  {
    # Brand new: create the family tree and add it a new classification:
    classif <- languageclassification();
    # Add a whole path including at least one internal node:
    newick.path <- rep("(",path.len);
    for( i in path.len:1 )
    {
      newick.path <- c(newick.path, path[i]);
      if( !is.null(brlens) ) newick.path <- c(newick.path, ":", brlens[i]);
      newick.path <- c(newick.path, ")");
    }
    newick.path <- c(newick.path, ";")
    new.tree <- familytree(text=paste0(newick.path,collapse=""), tree.name=path[1]);
    classif[[length(classif)+1]] <- new.tree;
    # Return it:
    return (classif);
  } else if( !inherits(classif, "languageclassification") )
  {
    return (NULL);
  } else
  {
    # Try to find the language family (if it already exists in the classification):
    found <- FALSE;
    for( i in 1:length(classif) )
    {
      if( get.name(classif[[i]]) == path[1] )
      {
        # Found!
        found <- TRUE;
        # Add the path the existing tree:
        new.tree <- add.tree.path(classif[[i]], path, brlens, warn.on.duplicated.tip);
        # and replace the tree:
        if( !is.null(new.tree) ) classif[[i]] <- new.tree;
        break;
      }
    }
    if( !found )
    {
      # Create a new family and add it to the classification:
      # Add a whole path including at least one internal node:
      newick.path <- rep("(",path.len);
      for( i in path.len:1 )
      {
        newick.path <- c(newick.path, path[i]);
        if( !is.null(brlens) ) newick.path <- c(newick.path, ":", brlens[i]);
        newick.path <- c(newick.path, ")");
      }
      newick.path <- c(newick.path, ";")
      new.tree <- familytree(text=paste0(newick.path,collapse=""), tree.name=path[1]);
      if( !is.null(new.tree) ) classif[[length(classif)+1]] <- new.tree;
    }
    # Return the classification:
    return (classif);
  }
}
add.language.to.languageclassification <- cmpfun(add.language.to.languageclassification);
# classif <- NULL;
# classif <- add.language.to.languageclassification(classif, c("f1","a","b"));
# classif <- add.language.to.languageclassification(classif, c("f2","c","d","e"));
# classif <- add.language.to.languageclassification(classif, c("f1","a","f","g"));
# classif <- add.language.to.languageclassification(classif, c("f1","a","x"));
# classif <- add.language.to.languageclassification(classif, c("f1","a","g")); # adding duplicated language

# Export a set of trees to a table containing the tree names and the Newick format:
export.languageclassification <- function( roots,                      # the language families
                                           dir.name="./",              # the destination directory
                                           classification="",          # the classification
                                           export.nexus=TRUE,          # export in the Nexus format?
                                           nexus.translate.block=TRUE, # use a translate block in the Nexus format?
                                           export.csv=TRUE,            # export in CSV (TAB-separated) format?
                                           method="none",              # the method
                                           constant=1.0,               # the constant (for those methods that need one)
                                           distmatrix=NA,              # the distance matrix (for those methods that need one)
                                           replace.NA.brlen.with=NA,   # some methods produce NA branch length (if the actual brlen cannot be estimated from the data): replace it by another numeric value?
                                           restore.collapsed.singles=TRUE, # some standard methods need to collapse singles, but we can restore them
                                           filename.postfix="",        # the suffix to the attached to the results file name
                                           quotes="'",                 # quotes
                                           parallel.mc.cores=1        # multicore processing?
                                         )
{
  if( !inherits(roots, "languageclassification") )
  {
    cat( "Can only export an object of class FamilyTrees!\n" );
    return (FALSE);
  }
  
  # The filename:
  if( classification != "" )
  {
    folder.name <- paste( dir.name, "/", classification, "/", sep="" );
  } else
  {
    folder.name <- dir.name;
  }
  if( !file.exists( folder.name ) ) dir.create( folder.name );
  if( classification != "" )
  {
    filename <- paste( folder.name, classification, "-newick", filename.postfix, ".csv", sep="" );
    filename.nex <- paste( folder.name, classification, "-nexus", filename.postfix, ".nex", sep="" );
  } else
  {
    filename <- paste( folder.name, "newick", filename.postfix, ".csv", sep="" );
    filename.nex <- paste( folder.name, "nexus", filename.postfix, ".nex", sep="" );
  }
  
  if( method == "none" )
  {
    # No brlengh: simply export the trees:

    # The translate block for NEXUS files:
    translate.block <- NULL; 
    if( export.nexus && nexus.translate.block )
    {
      translate.block <- lapply( 1:length(roots), function(i){
            root <- roots[[i]];
            if( !is.null(root) )
            {
              nodes.terminal <- extract.languages(root);
              nodes.internal <- extract.internal.nodes(root);
              ret.val <- data.frame( "Name"=c( nodes.terminal, nodes.internal ), "Internal"=c( rep(FALSE,length(nodes.terminal)), rep(TRUE,length(nodes.internal)) ), "ID"=NA );
              return (ret.val);
            }
        } );
      translate.block <- do.call(rbind, translate.block); translate.block <- unique(translate.block); translate.block$Name <- as.character(translate.block$Name);
      translate.block$Name[ translate.block$Name == "" ] <- "''"; # make sure the empty nodes appear as such
      translate.block$ID <- paste("n",1:nrow(translate.block),sep="");
    }

    # CSV file:
    if( export.csv )
    {
      cat( "Family\tTree\n", file=filename, append=FALSE );
    }
    # NEXUS file:
    if( export.nexus )
    {
      cat( "#NEXUS\n", file=filename.nex, append=FALSE );
      cat( "\n[Automatically generated by export.languageclassification() in FamilyTrees.R, (c) Dan Dediu 2014-2015.]\n\n", file=filename.nex, append=TRUE );
      cat( "begin trees;\n", file=filename.nex, append=TRUE );
      if( nexus.translate.block && !is.null(translate.block) )
      {
        # Export the translate block:
        cat( "\ttranslate\n", file=filename.nex, append=TRUE );
        for( i in 1:nrow(translate.block) )
        {
          cat( paste( "\t\t", translate.block$ID[i], "\t", translate.block$Name[i], ifelse(i < nrow(translate.block), ",", ";"), "\n", sep="" ), file=filename.nex, append=TRUE );
        }
        cat( "\n", file=filename.nex, append=TRUE );
      }
    }

    # The trees:
    for( i in 1:length(roots) )
    {
      root <- roots[[i]];
      if( !is.null(root) )
      {
        if( export.csv )
        {
          cat( paste( get.name(root), "\t", as.character( root ), "\n", sep="" ), file=filename, append=TRUE );
        }
        if( export.nexus )
        {
          s <- as.character( root );
          if( export.nexus && !is.null(translate.block) )
          {
            # Replace the names with their IDs:
            s.split <- strsplit(s, quotes, fixed=TRUE)[[1]];
            s.new <- NULL;
            for( k in 1:length(s.split) )
            {
              if( k %% 2 == 1 )
              {
                s.new <- c(s.new, s.split[k] );
              } else
              {
                s.new <- c(s.new, translate.block$ID[ translate.block$Name == paste( quotes, s.split[k], quotes, sep="" ) ] );
              }
            }
            s <- paste( s.new, collapse="", sep="" );
          }
          cat( paste( "\ttree ", nexus.string(get.name(root), force.to.nexus=TRUE), " = ", s, "\n", sep="" ), file=filename.nex, append=TRUE );
        }
      }
    }   

    if( export.nexus )
    {
      cat( "end;\n\n", file=filename.nex, append=TRUE );
    }
  } else
  {
    # Apply the brlength method and export the trees:
    tmp <- mclapply( 1:length(roots), function(i) 
                    { 
                      if( PRINT_DEBUG_MESSAGES ) cat( "    Family ", get.name(roots[[i]]), "(", i, " out of ", length(roots), "): " );
                      tree <- compute.brlen( roots[[i]], 
                                             method=method, 
                                             constant=constant, 
                                             distmatrix=distmatrix,
                                             replace.NA.brlen.with=replace.NA.brlen.with,
                                             restore.collapsed.singles=restore.collapsed.singles); 
                      if( !is.null(tree$tree) )
                      {
                        if( PRINT_DEBUG_MESSAGES ) cat( "OK\n" );
                        return ( list( "name"=get.name(roots[[i]]), "tree"=tree$tree, "comment"=tree$comment, "success"=TRUE ) );
                      } else
                      {
                        #stop("ERROR\n");
                        return ( list( "name"=get.name(roots[[i]]), "tree"=NULL, "comment"=tree$comment, "success"=FALSE ) );
                      }
                    }, mc.cores=parallel.mc.cores
                  );
  
    # The translate block for NEXUS files:
    translate.block <- NULL; 
    if( export.nexus && nexus.translate.block )
    {
      translate.block <- lapply( 1:length(tmp), function(i){            
            if( !is.null(tmp[[i]]) && (class(tmp[[i]]) == "list") && 
                  ("success" %in% names(tmp[[i]])) && ("tree" %in% names(tmp[[i]])) && 
                  !is.null(tmp[[i]]$success) && (tmp[[i]]$success == TRUE) && !is.null(tmp[[i]]$tree) )
            {
              root <- tmp[[i]]$tree;
              nodes.terminal <- extract.languages(root);
              nodes.internal <- extract.internal.nodes(root);
              ret.val <- data.frame( "Name"=c( nodes.terminal, nodes.internal ), "Internal"=c( rep(FALSE,length(nodes.terminal)), rep(TRUE,length(nodes.internal)) ), "ID"=NA );
              return (ret.val);
            } else
            {
              return (NULL);
            }
        } );
      translate.block <- do.call(rbind, translate.block); translate.block <- unique(translate.block); translate.block$Name <- as.character(translate.block$Name);
      if( !is.null(translate.block) )
      {
        translate.block$Name[ translate.block$Name == "" ] <- "''"; # make sure the empty nodes appear as such
        translate.block$ID <- paste("n",1:nrow(translate.block),sep="");
      }
    }

    # CSV file:
    if( export.csv )
    {
      cat( "Family\tSuccess\tComments\tTree\n", file=filename, append=FALSE );
    }

    # NEXUS file:
    if( export.nexus )
    {
      cat( "#NEXUS\n", file=filename.nex, append=FALSE );
      cat( "\n[Automatically generated by export.languageclassification() in FamilyTrees.R, (c) Dan Dediu 2014-2015.]\n\n", file=filename.nex, append=TRUE );
      cat( "begin trees;\n", file=filename.nex, append=TRUE );
      if( nexus.translate.block && !is.null(translate.block) )
      {
        # Export the translate block:
        cat( "\ttranslate\n", file=filename.nex, append=TRUE );
        for( i in 1:nrow(translate.block) )
        {
          cat( paste( "\t\t", translate.block$ID[i], "\t", translate.block$Name[i], ifelse(i < nrow(translate.block), ",", ";"), "\n", sep="" ), file=filename.nex, append=TRUE );
        }
        cat( "\n", file=filename.nex, append=TRUE );
      }
    }

    # The trees:
    for( i in 1:length(tmp) )
    {
      if( !is.null(tmp[[i]]) && (class(tmp[[i]]) == "list") && 
                  ("success" %in% names(tmp[[i]])) && ("tree" %in% names(tmp[[i]])) &&  ("name" %in% names(tmp[[i]])) &&   ("comment" %in% names(tmp[[i]])) && 
                  !is.null(tmp[[i]]$success) )
      {
        if( export.csv )
        {
          cat( paste( tmp[[i]]$name, "\t", ifelse(tmp[[i]]$success==TRUE,"SUCCESS","FAIL"), "\t", tmp[[i]]$comment, "\t", tmp[[i]]$tree, "\n", sep="" ), file=filename, append=TRUE );
        }
        if( export.nexus && (tmp[[i]]$success == TRUE) && !is.null(tmp[[i]]$tree) )
        {
          s <- as.character( tmp[[i]]$tree );
          if( export.nexus && !is.null(translate.block) )
          {
            # Replace the names with their IDs:
            s.split <- strsplit(s, quotes, fixed=TRUE)[[1]];
            s.new <- NULL;
            for( k in 1:length(s.split) )
            {
              if( k %% 2 == 1 )
              {
                s.new <- c(s.new, s.split[k] );
              } else
              {
                s.new <- c(s.new, translate.block$ID[ translate.block$Name == paste( quotes, s.split[k], quotes, sep="" ) ] );
              }
            }
            s <- paste( s.new, collapse="", sep="" );
          }
          cat( paste( "\ttree ", nexus.string(tmp[[i]]$name, force.to.nexus=TRUE), " = ", s, "\n", sep="" ), file=filename.nex, append=TRUE );
        }
      }
    }

    if( export.nexus )
    {
      cat( "end;\n\n", file=filename.nex, append=TRUE );
    }
  }
  
  return (TRUE);
}
export.languageclassification <- cmpfun(export.languageclassification);


