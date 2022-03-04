deme.decomp <- function(ED, n.deme, node.indices){
  event.times <- sort(unique(ED[,6]))
  time.increments <- diff(event.times)

  k <- matrix(0, nrow = length(event.times) - 1, ncol = n.deme)
  root.row <- node.indices[ED[is.na(ED[,2]),1]]
  active.nodes <- ED[root.row, 3:4]
  k[1,ED[root.row, 5]] <- 2

  for (i in 2 : (length(event.times) - 1)){
    k[i,] <- k[i-1,]
    active.rows <- node.indices[active.nodes]
    current.indices <- which(ED[active.rows, 6] == event.times[i])
    current.rows <- active.rows[current.indices]

    if (length(current.rows) > 1){ #Multiple leaves
      for (j in current.rows){
        current.deme <- ED[j, 5]
        k[i, current.deme] <- k[i, current.deme] - 1
      }
    } else{
      if (anyNA(ED[current.rows, 3])){ #Leaf
        current.deme <- ED[current.rows, 5]
        k[i, current.deme] <- k[i, current.deme] - 1
      } else if (!anyNA(ED[current.rows, 4])){ #Coalescence
        current.deme <- ED[current.rows, 5]
        k[i, current.deme] <- k[i, current.deme] + 1
        active.nodes <- c(active.nodes, ED[current.rows, 3:4])
      } else{ #Migration
        current.deme <- ED[current.rows, 5]
        k[i, current.deme] <- k[i, current.deme] - 1
        current.child <- node.indices[ED[current.rows, 3]]
        child.deme <- ED[current.child, 5]
        k[i, child.deme] <- k[i, child.deme] + 1
        active.nodes <- c(active.nodes, ED[current.rows, 3])
      }
    }
    active.nodes <- active.nodes[-current.indices]
  }

  return(k)
}