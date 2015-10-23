#################################################################################################################################
#
# R script for computing a generalization of the "genetic distance" between languages as defined originally in 
#    Maurits, L. & Griffiths, T.L. (2015) Tracing the roots of syntax with Bayesian phylogenetics, PNAS 111 (37):13576--13581
# section Materials and Methods, subsection Constructing trees, page 13580.
# The generalization is to allow any hierarchical classification (not just the Ethnologue) and 
# any alpha (by default set to the reportoed value of 0.69); for languages from diferent families the distance is set to NA.
#
# The method is described in the paper as:
#
# The genetic method uses genetic classifications of languages taken from Ethnologue. 
# The classifications are used to assign a "genetic distance" between languages as follows. 
# Let language A be classified as A1 ⊃ A2 ⊃ . . . An and language B as B1 ⊃ B2 ⊃ . . . Bm. 
# Without loss of generality, let n ≤ m.
# We discarded Bn+1, . . ., Bm if necessary. 
# If Ai = Bi for i = 1, . . ., k, i.e., A and B are classified identically for the first k refinements of their family, then the distance is 
#     M − sum(i=1..k, α^i), 
# where α = 0.69 and M is the maximum value the sum can take, i.e., when k = n, so that identically classified languages have a distance of zero. 
# Since α < 1, each additional matching refinement contributes less to the sum, reflecting the fact that later refinements tend to be much more fine grained. 
#
#
# The description was later refiend in an e-mail exchange (2 June 2015) with Luke Maurits as follows:
#
# If I've understood your example tree correctly, then I do think you have followed the Genetic method correctly.
# To answer your first question, it's not actually the case that A = A1 or A = An.
# In this notation, A represents the language (say A = "Lithuanian"), A1 is the highest-level classification in Ethnologue (e.g. A1 = "Indo-European") 
# and An is the lowest-level classification (e.g. A4 = "Eastern", which should be interpreted as "Eastern Baltic" as A3 = "Baltic").
# But I can see that this system is not at all clear from the paper, which I
# apologise for.  You were spot on with the space pressures!
# 
# [The example I used in my request for explanation was:
#   root - S1 -S2 - S3.1 - S4.1 - A
#                   +-  S3.2 - B
#                   +-  S3.3 - C
# ]
#   
# This is equivalent to "discarding the leaves" as you asked about in your example, so the different "paths" would be:
# 
# A: root, S1, S2, S3.1, S4.1
# B: root, S1, S2, S3.2
# C: root, S1, S2, S3.3
# 
# As mentioned ("We discarded Bn+1, . . ., Bm if necessary"), these paths are truncated to the shortest of the pair when computing a distance, so we would have:
# 
# A: root, S1, S2, S3.1
# B: root, S1, S2, S3.2
# C: root, S1, S2, S3.3
# 
# So we do indeed find that d(A,B) = d(A,C) = d(B,C), i.e. all three languages are equidistant from one another.
# I agree that this might seem unexpected, since the full classification plainly tells us that B and C are closer to eachother than either is to A. 
# However, remembering that the purpose of each of the three methods is simply to generate a "reasonable" tree to use as a kind of central point for a
# region of treespace to sample from, truncating the classifications makes the method much simpler and should not (I think!) introduce anything too anomolous.
# Notice that if language A had a sibling language D which was also descended from S4.1, then A and D would have a distance from each other of zero, 
# which is less than the distances betwen A and B, A and C, D and A and D and B, so there should still end up being some internal structure in the (A,B,C,D) clade.
#
#################################################################################################################################

# Load the tree processing helper functions:
source("../../../code/FamilyTrees.R");

# Compute the MG2015 distances between all pairs of languages given a classification (as a FamilyTrees object) and alpha value
MG2015.distances <- function( classif, alpha=0.69, interfamily.distance=NA )
{
  # Preconditions:
  if( !inherits(classif,"languageclassification") || length(classif) == 0 )
  {
    cat("The langauge classification is of the wrong type or empty!\n");
    return (NULL);
  }
  
  # Build the list of all languages in the classification:
  lgs <- unlist(lapply(classif, extract.languages));
  if( length(lgs) < 1 )
  {
    cat("There must be at least one language!\n");
    return (NULL);    
  }
  # Initialize the distances matrix:
  m <- matrix(interfamily.distance, nrow=length(lgs), ncol=length(lgs)); rownames(m) <- colnames(m) <- lgs;
  
  # For a given classification and two language names, compute the MG2015 distance between the languages using the given alpha:
  MG2015.dist <- function(fam, lg1, lg2, alpha)
  {
    if( lg1 == lg2 ) return (0); # distance between language and self is 0
    # p1 and p2 are the paths from the root to the languages
    p1 <- extract.path(fam, lg1); p2 <- extract.path(fam, lg2);
    if( is.null(p1) || is.null(p2) ) return (NA); # distance is only defined if both languages belong to this tree
    # the lengths of the two paths:
    l1 <- length(p1); l2 <- length(p2);
    if( l1 < 1 || l2 < 1 ) return (NA); # distance is only defined if both languages belong to this tree
    # the potentially shared path length:
    n <- min(l1, l2);
    # the maximum possible sum for such a path:
    M <- sum( vapply( 1:n, function(i) alpha^i, numeric(1) ) );
    # the length k of the shared path between the two languages:
    k <- which.min( vapply( 1:n, function(i) p1[i] == p2[i], logical(1) ) ) - 1;
    # and the distance:
    d <- M - ifelse( k <= 0, 0, sum( vapply( 1:k, function(i) alpha^i, numeric(1) ) ) );
    # return it:
    return (d);
  }
  
  # Compute the distances only between languages from the same classification (between classifications the distance is by definition interfamily.distance):
  cat(paste0("MG2015 distance for classification ", get.name(classif), " (", length(classif), " families):\n"));
  for(i in seq_along(classif)) # i = grep("Indo",get.family.names(classif),fixed=TRUE)
  {
    fam <- classif[[i]]; # the i-th family
    cat(paste0("   family ", get.name(fam), " (", i, " of ", length(classif), ")\n"));
    if( !inherits(fam,"familytree") || is.null(fam) ) next;
    fam.lgs <- extract.languages(fam);
    if( length(fam.lgs) < 1 ) next;
    m[fam.lgs, fam.lgs] <- 0; # initialize the within-family distances to 0
    if( length(fam.lgs) < 2 ) next;
    for(j in 1:(length(fam.lgs)-1)) # j = grep("[i-eng]",fam.lgs,fixed=TRUE)
    {
      for(k in (j+1):length(fam.lgs)) # k = grep("[i-deu]",fam.lgs,fixed=TRUE)
      {
        d <- MG2015.dist(fam=fam, lg1=fam.lgs[j], lg2=fam.lgs[k], alpha=alpha); # not the optimum way to compute, but very clear
        m[fam.lgs[j], fam.lgs[k]] <- m[fam.lgs[k], fam.lgs[j]] <- d;
      }
    }
  }
  
  # Return the distances matrix:
  return (m);
}
MG2015.distances <- cmpfun(MG2015.distances);
## TESTS:
# m <- MG2015.distances(classif=languageclassification(classif.name="Ethnologue", csv.file="../../../output/ethnologue/ethnologue-newick.csv", csv.name.column="Family", csv.tree.column="Tree"), alpha=0.69, interfamily.distance=NA)


# Compute the distances for all four classifications:
alpha <- 0.69;

# WALS:
MG2015.WALS <- MG2015.distances(classif=languageclassification(classif.name="WALS", 
                                                               csv.file=paste0("../../../output/wals/wals-newick.csv"), 
                                                               csv.name.column="Family", csv.tree.column="Tree"), 
                                alpha=alpha, interfamily.distance=NA);
save(MG2015.WALS, file=paste0("./","MG2015-wals-alpha=",alpha,".RData"), compress="xz");

# Ethnologue:
MG2015.Ethnologue <- MG2015.distances(classif=languageclassification(classif.name="Ethnologue", 
                                                               csv.file=paste0("../../../output/ethnologue/ethnologue-newick.csv"), 
                                                               csv.name.column="Family", csv.tree.column="Tree"), 
                                alpha=alpha, interfamily.distance=NA);
save(MG2015.Ethnologue, file=paste0("./","MG2015-ethnologue-alpha=",alpha,".RData"), compress="xz");

# Glottolog:
MG2015.Glottolog <- MG2015.distances(classif=languageclassification(classif.name="Glottolog", 
                                                                     csv.file=paste0("../../../output/glottolog/glottolog-newick.csv"), 
                                                                     csv.name.column="Family", csv.tree.column="Tree"), 
                                      alpha=alpha, interfamily.distance=NA);
save(MG2015.Glottolog, file=paste0("./","MG2015-glottolog-alpha=",alpha,".RData"), compress="xz");

# AUTOTYP:
MG2015.AUTOTYP <- MG2015.distances(classif=languageclassification(classif.name="AUTOTYP", 
                                                                    csv.file=paste0("../../../output/autotyp/autotyp-newick.csv"), 
                                                                    csv.name.column="Family", csv.tree.column="Tree"), 
                                     alpha=alpha, interfamily.distance=NA);
save(MG2015.AUTOTYP, file=paste0("./","MG2015-autotyp-alpha=",alpha,".RData"), compress="xz");
