#####################################################################################
# GetBipEnergy is the function to calculate the energy of a given bipartite netowrk	#
# which is defined as -\sum_{ij}\sum_{Neighbor(i'j')}(2*x_{ij}-1)*(2*x_{i'j'}-1)	#
# Mirroring needed on the lattice edges												#
# Input --> bipartite matrix, 0/1; Output --> energy, numeric						#
#####################################################################################

GetBipEnergy <- function(x){			# Bipartite matrix x
	
	n1 <- nrow(x)						# The number of the row-subject and col-subject
	n2 <- ncol(x)
	
	# Define a new matrix x.new which mirrors the elements at lattice edges
		
	x.new <- matrix(0, nrow = (n1+2), ncol = (n2+2))
	x.new[2:(n1+1), 2:(n2+1)] <- x		# The inner elements are mirrored from the given x
	x.new[1, 2:(n2+1)] <- x[2,]			# All the peripheral elements being flipped
	x.new[2:(n1+1), 1] <- x[,2]
	x.new[(n1+2), 2:(n2+1)] <- x[(n1-1),]
	x.new[2:(n1+1), (n2+2)] <- x[,(n2-1)]
	
	eng.mat <- matrix(0, nrow = n1, ncol = n2)
										# A matrix saving the energy contributions
	
	for(i in 2:(n1+1)){
		for (j in 2:(n2+1)){			# Calculate only for the half of the elements
			
			top <- x.new[(i-1), j]		# Find the four neighbors of (i,j) entry
			bottom <- x.new[(i+1), j]
			left <- x.new[i, (j-1)]
			right <- x.new[i, (j+1)]
			
			neighbor <- c(top, bottom, left, right)
			eng.mat[(i-1), (j-1)] <- -sum((2*x.new[i,j] - 1)*(2*neighbor - 1))
										# The energy contribution from (i,j) entry
		}
	}
	
	energy <- sum(eng.mat)				# Energy is the of enrgy contribution
	return(energy)						# Output the energy	
}
# End of GetBipEnergy

#####################################################################################
# GetEnergy.byRow is the function to calculate the energy contribution of a given	# 
# row from a binary matrix. Mirroring needed on the lattice edges					#
# Input: x --> bipartite matrix, binary; row --> specific rows, numeric				#
# Output --> energy contribution of the rows, numeric								#
#####################################################################################

GetEnergy.byRow <- function(x, row){
	
	n1 <- nrow(x)						# The number of the row-subject and col-subject
	n2 <- ncol(x)	
	n.row <- length(row)				# The number of rows to be calculated
	
	
	# Define the new matrix x.new
		
	x.new <- matrix(0, nrow = (n1+2), ncol = (n2+2))
	x.new[2:(n1+1), 2:(n2+1)] <- x		# The inner elements are mirrored from the given x
	x.new[1, 2:(n2+1)] <- x[2,]			# All the peripheral elements being flipped
	x.new[2:(n1+1), 1] <- x[,2]
	x.new[(n1+2), 2:(n2+1)] <- x[(n1-1),]
	x.new[2:(n1+1), (n2+2)] <- x[,(n2-1)]
	
	energy.row <- matrix(0, nrow = n.row, ncol = n2)
	for(r in 1:n.row){
		for(k in 2:(n2+1)){
			top <- x.new[row[r], k]			# Find the four neighbors of (i,j) entry
			bottom <- x.new[(row[r]+2), k]
			left <- x.new[(row[r]+1), (k-1)]
			right <- x.new[(row[r]+1), (k+1)]
			
			neighbor <- c(top, bottom, left, right)
			energy.row[r, k-1] <- -sum((2*x.new[(row[r]+1),k] - 1)*(2*neighbor - 1))
		}
	}
	
	return(rowSums(energy.row))
}
# End of GetEnergy.byRow

#####################################################################################
# EnergyUpdate is the function to calculate the energy change by permuting a pair of# 
# rows in a binary matrix. 															#
# Input: x --> matrix, binary; row1/row2 --> specific rows, numeric					#
# Output --> energy change, numeric													#
#####################################################################################

EnergyDiff <- function(x, row1, row2){
	
	n <- nrow(x)						# Nrow of x
	row.idx <- unique(c((row1-1):(row1+1), (row2-1):(row2+1)))
	row.idx <- row.idx[row.idx>=1 & row.idx <=n]
										# Rows that change the energy from permutation
	pmt <- 1:n
	pmt[row1] <- row2
	pmt[row2] <- row1
	x.pmt <- x[pmt,]					# matrix with permuted rows
	
	diff <- sum(GetEnergy.byRow(x.pmt, row.idx)) - sum(GetEnergy.byRow(x, row.idx))
	return(diff)
}
#End of EnergyDiff

#########################################################################################
# SA1 is the function permuting the row-nodes in row-blocks for searching ground state	#
# given a FIXED hierarchical tree structure												#
# Input: x --> Bipartite matrix, binary;												#
# Input: IDlist --> List of clusters containing ids in each of them,list;				#
# Input: ord --> Initial order of the tree nodes, numeric;								#
# Input: Ite -- Number of iterations, numeric.											#
# Output: A list with permuted x and updated order.										#
#########################################################################################

SA1 <- function(x, IDlist, ord, Ite = 1000){
	
	n <- nrow(x)										# Number of row-nodes
	nclust <- length(IDlist)							# Number of clusters
	clust.size <- unlist(lapply(IDlist, length))		# Cluster sizes
	
	p <- (clust.size-1)^2/sum((clust.size-1)^2)
				# Probability vector of permuting any pair in the cluster
				# Probability proportional to the cluster size square
				# Singletons never get picked (no permutations needed) 
	
	energy.best <- GetBipEnergy(x)
	# record <- energy.best
	
	for(j in 1:Ite){
		
		k <- sample(1:nclust, size = 1, prob = p)       # Sample one cluster
		node <- sample(IDlist[[k]], size = 2)			# Sample two nodes in the cluster
		idx1 <- which(ord == node[1])					# Locate the two nodes
		idx2 <- which(ord == node[2])
		
		pmt <- 1:n										# Reset the permuting pair
		pmt[idx1] <- idx2								# Permute the pair
		pmt[idx2] <- idx1
		x.new <- x[pmt,]								# Neighbor state to which it moves
		
		update <- EnergyDiff(x, idx1, idx2)				# Energy change calculated
		if(update < 0){
			x <- x.new									# Update the matrix
			ord <- ord[pmt]								# Update the order
			energy.best <- energy.best + update			# Update the best energy
			print(update)
		}
	}
		
	return(list(x,ord,energy.best))
}
# End of SA1

#################################################################################
# Optimize is the function searching for GLOBAL optimal energy given a list of	#
# permissible hierarchical tree structures; For each structure, SA1 is further	#
# applied to search for a LOCAL optimum; All local minimum is recorded			#									
# Input: x --> Bipartite matrix, binary;										#
# Input: idList --> List of clusters,list; Passed to SA1						#
# Input: clustStr --> List of permissible cluster-order, which conforms to the	#
# tree structure returned from DCG												#
# Output: A list with permuted x, optimal order & optimal energy.				#
#################################################################################

Optimize <- function(x, idList, clustStr){
	
	nclust <- length(idList)			# Number of clusters
	npermn <- length(clustStr)			# Number of permissible cluster permutation
	
	matrix.opt <- list()				# Set the optimum list
	for(m in 1:npermn){
		layout <- clustStr[[m]]
		order <- numeric(0)
		
		for(j in 1:length(layout)){
			order <- c(order, idList[[layout[j]]])
		}								# Initial row-nodes order
		
		x.int <- x[order,]				# Initial matrix
		output <- SA1(x.int, idList, order, Ite = 2000)
										# Searching for the LOCAL optimum
		matrix.opt[[m]] <- list(x.int, output[[1]], output[[2]], output[[3]])
	}
	
	return(matrix.opt)
}
# End of SA2 
