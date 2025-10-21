


# -------------------------- Voting rules available --------------------------
rule_map <- c(
  borda          = "Borda",
  plurality      = "Plurality",
  approval       = "Approval (≥ mean)",
  approval2      = "Approval (≥ 5)",
  mean_winner    = "Mean Utility",
  copeland       = "Copeland",
  general_borda  = "Gen. Borda (top-q)",
  general_borda_3  = "Gen. Borda (top 3)",
  general_borda_5  = "Gen. Borda (top 5)",
  multiple_votes   = "k-Approval (top-g)",
  multiple_votes_3 = "k-Approval (top 3)",
  multiple_votes_5 = "k-Approval (top 5)",
  ostracism        = "Ostracism (+1/-1)",
  ostracism_3      = "Ostracism (↑3, ↓3)",
  ostr_gen_borda   = "Ostra-Borda (±3..±1)",
  condorcet      = "Condorcet",
  dictator       = "Dictator",
  random         = "Random",
  majority       = "Majority (>50% first prefs)",
  two_round      = "Two-Round (runoff)",
  irv            = "Instant Runoff (RCV)",
  typical_judgment  = "Typical Judgment",
  majority_judgment = "Majority Judgment",
  dictator_r          = "Dictator (ranked)",
  random_r            = "Random (ranked)",
  majority_r          = "Majority (ranked)",
  two_round_r         = "Two-Round (ranked)",
  irv_r               = "Instant Runoff (ranked)",
  condorcet_r         = "Condorcet (Copeland ranking)",
  typical_judgment_r  = "Typical Judgment (ranked)",
  majority_judgment_r = "Majority Judgment (ranked)"
)


# ---------------------- Voting rules implemented -------------------------
voting_strategy <- function(utilities, scheme = "borda", q = NULL, g = NULL, dictator = NULL) {
  utilities <- as.matrix(utilities)
  n_voters <- nrow(utilities); n_alternatives <- ncol(utilities)
  alt_names <- colnames(utilities)
  
  make_out <- function(scores, ballots = NULL, winner_idx = NULL, picked = NULL) {
    names(scores) <- alt_names
    list(scores = as.numeric(scores),
         ballots = ballots,
         winner_set = if (is.null(winner_idx)) which(scores == max(scores, na.rm=TRUE)) else winner_idx,
         picked = if (is.null(picked)) sample(which(scores == max(scores, na.rm=TRUE)), 1) else picked)
  }
  
  # default empty ballots
  ballots <- matrix(0, nrow = n_voters, ncol = n_alternatives,
                    dimnames = list(NULL, alt_names))
  
  if (scheme == "borda") {
    ballots <- t(apply(utilities, 1, function(x) {
      r <- rank(-x, ties.method = "random"); n_alternatives - r
    }))
    return(make_out(colSums(ballots), ballots))
    
  } else if (scheme == "plurality") {
    ballots <- t(apply(utilities, 1, function(x) {
      ties <- which(x == max(x)); v <- integer(length(x)); v[sample(ties,1)] <- 1L; v
    }))
    return(make_out(colSums(ballots), ballots))
    
  } else if (scheme == "approval") {
    ballots <- t(apply(utilities, 1, function(x) {
      v <- as.integer(x >= mean(x, na.rm = TRUE)); if (!sum(v)) v[which.max(x)] <- 1L; v
    }))
    return(make_out(colSums(ballots), ballots))
    
  } else if (scheme == "approval2") {
    ballots <- t(apply(utilities, 1, function(x) {
      v <- as.integer(!is.na(x) & x >= 5); if (!sum(v)) v[which.max(x)] <- 1L; v
    }))
    return(make_out(colSums(ballots), ballots))
    
  } else if (scheme == "mean_winner") {
    s <- colMeans(utilities, na.rm = TRUE)
    return(make_out(s, ballots = NULL))
    
  } else if (scheme == "copeland") {
    wins <- matrix(0L, n_alternatives, n_alternatives)
    for (i in 1:(n_alternatives-1)) for (j in (i+1):n_alternatives) {
      ij <- sum(utilities[,i] > utilities[,j], na.rm = TRUE)
      ji <- sum(utilities[,j] > utilities[,i], na.rm = TRUE)
      if (ij > ji) {
        wins[i, j] <- 1
      } else if (ji > ij) {
        wins[j, i] <- 1
      } else {
        wins[i, j] <- 0.5
        wins[j, i] <- 0.5
      }
    }
    s <- rowSums(wins)
    return(make_out(s, ballots = NULL))
    
  } else if (scheme == "general_borda") {
    if (is.null(q) || !is.finite(q)) q <- n_alternatives
    q <- max(1L, min(as.integer(q), n_alternatives))
    ballots <- t(apply(utilities, 1, function(x) {
      r <- rank(-x, ties.method = "random"); pts <- integer(n_alternatives)
      top <- which(r <= q); pts[top] <- q - r[top] + 1L; pts
    }))
    return(make_out(colSums(ballots), ballots))
    
  } else if (scheme == "general_borda_3") {
    q <- 3L
    q <- max(1L, min(q, n_alternatives))
    ballots <- t(apply(utilities, 1, function(x) {
      r <- rank(-x, ties.method = "random"); pts <- integer(n_alternatives)
      top <- which(r <= q); pts[top] <- q - r[top] + 1L; pts
    }))
    return(make_out(colSums(ballots), ballots))
    
  } else if (scheme == "general_borda_5") {
    q <- 5L
    q <- max(1L, min(q, n_alternatives))
    ballots <- t(apply(utilities, 1, function(x) {
      r <- rank(-x, ties.method = "random"); pts <- integer(n_alternatives)
      top <- which(r <= q); pts[top] <- q - r[top] + 1L; pts
    }))
    return(make_out(colSums(ballots), ballots))
    
  } else if (scheme == "multiple_votes") {
    if (is.null(g) || !is.finite(g)) g <- 1L
    g <- max(1L, min(as.integer(g), n_alternatives))
    ballots <- t(apply(utilities, 1, function(x) {
      r <- rank(-x, ties.method = "random"); v <- integer(n_alternatives); v[r <= g] <- 1L; v
    }))
    return(make_out(colSums(ballots), ballots))
    
  } else if (scheme == "multiple_votes_3") {
    g <- 3L
    g <- max(1L, min(g, n_alternatives))
    ballots <- t(apply(utilities, 1, function(x) {
      r <- rank(-x, ties.method = "random"); v <- integer(n_alternatives); v[r <= g] <- 1L; v
    }))
    return(make_out(colSums(ballots), ballots))
    
  } else if (scheme == "multiple_votes_5") {
    g <- 5L
    g <- max(1L, min(g, n_alternatives))
    ballots <- t(apply(utilities, 1, function(x) {
      r <- rank(-x, ties.method = "random"); v <- integer(n_alternatives); v[r <= g] <- 1L; v
    }))
    return(make_out(colSums(ballots), ballots))
    
  } else if (scheme == "ostracism") {
    ballots <- t(apply(utilities, 1, function(x) {
      v <- integer(length(x)); v[which.max(x)] <-  1L; v[which.min(x)] <- -1L; v
    }))
    return(make_out(colSums(ballots), ballots))
    
  
    } else if (scheme == "ostracism_3") {
    ballots <- t(apply(utilities, 1, function(x) {
      n <- length(x); k <- min(3L, n)
      r <- rank(-x, ties.method = "random")     # 1 = best
      v <- integer(n)
      top    <- which(r <= k)
      bottom <- which(r > n - k)
      if (length(top))    v[top]    <- v[top]    + 1L
      if (length(bottom)) v[bottom] <- v[bottom] - 1L
      v
    }))
    return(make_out(colSums(ballots), ballots))
  
    } else if (scheme == "ostr_gen_borda") {
    ballots <- t(apply(utilities, 1, function(x) {
      n <- length(x); k <- min(3L, n)
      r <- rank(-x, ties.method = "random")     # 1 = best
      v <- integer(n)
      
      # Top k: points k, k-1, ..., 1
      top <- which(r <= k)
      if (length(top)) v[top] <- k - r[top] + 1L
      
      # Bottom k: points -k, -(k-1), ..., -1
      bottom <- which(r > n - k)
      if (length(bottom)) v[bottom] <- v[bottom] - (r[bottom] - (n - k))
      
      v
    }))
    return(make_out(colSums(ballots), ballots))
  
    } else if (scheme == "condorcet") {
    wins <- matrix(0L, n_alternatives, n_alternatives)
    for (i in 1:(n_alternatives - 1)) for (j in (i + 1):n_alternatives) {
      ij <- sum(utilities[, i] > utilities[, j], na.rm = TRUE)
      ji <- sum(utilities[, j] > utilities[, i], na.rm = TRUE)
      if (ij > ji) wins[i, j] <- 1L else if (ji > ij) wins[j, i] <- 1L
    }
    s <- rowSums(wins)
    cw <- which(s == (n_alternatives - 1L))              # strict Condorcet winners
    wset <- if (length(cw)) cw else which(s == max(s))    # expose maximal set
    picked <- sample(wset, 1)
    # score as a 1 for picked (pure winner) WITHOUT inflating by n_voters:
    score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
    return(make_out(score_vec, ballots = NULL, winner_idx = wset, picked = picked))
    
  } else if (scheme == "dictator") {
    d <- if (!is.null(dictator) && is.finite(dictator) && dictator >= 1 && dictator <= n_voters)
      as.integer(dictator) else sample.int(n_voters, 1)
    ties <- which(utilities[d, ] == max(utilities[d, ], na.rm = TRUE))
    picked <- sample(ties, 1)
    score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
    return(make_out(score_vec, ballots = NULL, winner_idx = ties, picked = picked))
    
  } else if (scheme == "random") {
    picked <- sample.int(n_alternatives, 1)
    score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
    return(make_out(score_vec, ballots = NULL, winner_idx = picked, picked = picked))
    
  } else if (scheme == "majority") {
    # First-preference counts (plurality) with random tie-breaking at voter level
    first_choice <- apply(utilities, 1, function(x) {
      ties <- which(x == max(x, na.rm = TRUE)); sample(ties, 1)
    })
    counts <- tabulate(first_choice, nbins = n_alternatives)
    
    # Absolute majority threshold
    if (max(counts) > n_voters / 2) {
      winners <- which(counts == max(counts))
    } else {
      # No absolute majority -> choose plurality leader (random if tie)
      winners <- which(counts == max(counts))
    }
    picked <- sample(winners, 1)
    
    score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
    return(make_out(score_vec, ballots = NULL, winner_idx = winners, picked = picked))
    
  } else if (scheme == "two_round") {
    # Round 1: plurality
    first_choice <- apply(utilities, 1, function(x) {
      ties <- which(x == max(x, na.rm = TRUE)); sample(ties, 1)
    })
    counts1 <- tabulate(first_choice, nbins = n_alternatives)
    
    # Winner in round 1?
    if (max(counts1) > n_voters / 2) {
      winners <- which(counts1 == max(counts1))
      picked  <- sample(winners, 1)
      score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
      return(make_out(score_vec, ballots = NULL, winner_idx = winners, picked = picked))
    }
    
    # No majority: pick top two for runoff (handle ties by randomizing among top counts)
    top_count <- max(counts1)
    top_set   <- which(counts1 == top_count)
    if (length(top_set) > 2) {
      finalists <- sample(top_set, 2)
    } else if (length(top_set) == 2) {
      finalists <- top_set
    } else {
      # Unique top; need the runner-up (ties allowed)
      second_count <- sort(counts1, decreasing = TRUE)[2]
      runners <- which(counts1 == second_count)
      finalists <- c(top_set, sample(runners, 1))
    }
    
    a <- finalists[1]; b <- finalists[2]
    
    # Round 2: pairwise majority using utilities between a and b
    ab_pref <- sum(utilities[, a] > utilities[, b], na.rm = TRUE)
    ba_pref <- sum(utilities[, b] > utilities[, a], na.rm = TRUE)
    
    if (ab_pref > ba_pref) {
      picked <- a
    } else if (ba_pref > ab_pref) {
      picked <- b
    } else {
      picked <- sample(c(a, b), 1)  # exact tie
    }
    
    score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
    return(make_out(score_vec, ballots = NULL, winner_idx = finalists, picked = picked))
    
  } else if (scheme == "irv") {
    # Build a strict ranking per voter (random tie-breaks within each voter)
    ranks <- t(apply(utilities, 1, function(x) {
      # order descending by utility; break ties by random permutation of tied positions
      o <- order(-x, runif(length(x)))
      o
    }))  # rows: voters, columns: rank positions (1 = top)
    
    active <- rep(TRUE, n_alternatives)
    picked <- NA_integer_
    
    repeat {
      # First preferences among active candidates
      first_prefs <- apply(ranks, 1, function(ord) {
        # find first active candidate in this voter's order
        for (idx in ord) if (active[idx]) return(idx)
        return(NA_integer_)  # should not happen if at least one active
      })
      first_prefs <- first_prefs[!is.na(first_prefs)]
      counts <- tabulate(first_prefs, nbins = n_alternatives)
      
      total_active_ballots <- sum(counts)
      if (total_active_ballots == 0) {
        # all exhausted (shouldn't occur here); pick randomly among active
        picked <- sample(which(active), 1); break
      }
      
      # Majority?
      if (max(counts) > total_active_ballots / 2) {
        winners <- which(counts == max(counts))
        picked  <- sample(winners, 1)
        break
      }
      
      # Eliminate the lowest count among active
      min_count <- min(counts[active])
      losers <- which(active & counts == min_count)
      # break ties at the bottom randomly, eliminate one at a time
      to_eliminate <- sample(losers, 1)
      active[to_eliminate] <- FALSE
      
      # If only one remains active, they win
      if (sum(active) == 1) { picked <- which(active); break }
    }
    
    score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
    return(make_out(score_vec, ballots = NULL, winner_idx = picked, picked = picked))
    
  } else if (scheme == "typical_judgment") {
    # Best median filter, then maximize Pb - Pw
    meds <- apply(utilities, 2, function(col) median(col, na.rm = TRUE))
    if (all(is.na(meds))) {
      picked <- sample.int(n_alternatives, 1)
      score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
      return(make_out(score_vec, ballots = NULL, winner_idx = picked, picked = picked))
    }
    best_med <- max(meds, na.rm = TRUE)
    cand_set <- which(meds == best_med)
    
    pb <- pw <- numeric(length(cand_set))
    for (k in seq_along(cand_set)) {
      j <- cand_set[k]
      x <- utilities[, j]
      x <- x[!is.na(x)]
      if (!length(x)) { pb[k] <- 0; pw[k] <- 0 } else {
        m <- stats::median(x)
        n <- length(x)
        pb[k] <- sum(x >  m) / n
        pw[k] <- sum(x <  m) / n
      }
    }
    val <- pb - pw
    winners <- cand_set[which(val == max(val, na.rm = TRUE))]
    picked  <- sample(winners, 1)
    
    score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
    return(make_out(score_vec, ballots = NULL, winner_idx = winners, picked = picked))
    
  } else if (scheme == "majority_judgment") {
    # Best median filter, then attach s = Pb if Pb>Pw, else s = -Pw; maximize s
    meds <- apply(utilities, 2, function(col) median(col, na.rm = TRUE))
    if (all(is.na(meds))) {
      picked <- sample.int(n_alternatives, 1)
      score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
      return(make_out(score_vec, ballots = NULL, winner_idx = picked, picked = picked))
    }
    best_med <- max(meds, na.rm = TRUE)
    cand_set <- which(meds == best_med)
    
    s <- numeric(length(cand_set))
    for (k in seq_along(cand_set)) {
      j <- cand_set[k]
      x <- utilities[, j]
      x <- x[!is.na(x)]
      if (!length(x)) { s[k] <- 0 } else {
        m <- stats::median(x)
        n <- length(x)
        pb <- sum(x > m) / n
        pw <- sum(x < m) / n
        if (pb > pw) s[k] <- pb else if (pw > pb) s[k] <- -pw else s[k] <- 0
      }
    }
    winners <- cand_set[which(s == max(s, na.rm = TRUE))]
    picked  <- sample(winners, 1)
    
    score_vec <- rep(0, n_alternatives); score_vec[picked] <- 1
    return(make_out(score_vec, ballots = NULL, winner_idx = winners, picked = picked))
  } else if (scheme == "dictator_r") {
    d <- if (!is.null(dictator) && is.finite(dictator) && dictator >= 1 && dictator <= n_voters)
      as.integer(dictator) else sample.int(n_voters, 1)
    scores <- as.numeric(utilities[d, ])
    # winner: dictator’s top choice (break ties randomly)
    wset <- which(scores == max(scores, na.rm = TRUE))
    picked <- sample(wset, 1)
    return(make_out(scores, ballots = NULL, winner_idx = wset, picked = picked))
    
  } else if (scheme == "random_r") {
    # purely random total order (no informational content; for comparability only)
    scores <- runif(n_alternatives)
    picked <- which.max(scores)
    return(make_out(scores, ballots = NULL, winner_idx = picked, picked = picked))
    
  } else if (scheme == "majority_r") {
    # rank by first-preference proportions
    first_choice <- apply(utilities, 1, function(x) {
      ties <- which(x == max(x, na.rm = TRUE)); sample(ties, 1)
    })
    counts <- tabulate(first_choice, nbins = n_alternatives)
    prop1  <- counts / sum(counts)
    scores <- prop1
    # majority winner if >50%, else plurality leader(s)
    if (max(counts) > n_voters / 2) wset <- which(counts == max(counts)) else wset <- which(counts == max(counts))
    picked <- sample(wset, 1)
    return(make_out(scores, ballots = NULL, winner_idx = wset, picked = picked))
    
  } else if (scheme == "two_round_r") {
    # Round 1: plurality
    first_choice <- apply(utilities, 1, function(x) {
      ties <- which(x == max(x, na.rm = TRUE)); sample(ties, 1)
    })
    counts1 <- tabulate(first_choice, nbins = n_alternatives)
    prop1   <- counts1 / sum(counts1)
    
    # Winner in round 1?
    if (max(counts1) > n_voters / 2) {
      winners <- which(counts1 == max(counts1))
      picked  <- sample(winners, 1)
      # ranking = promote this winner; others by prop1
      M <- 2.0  # larger than any prop1 gap
      scores <- prop1; scores[picked] <- max(prop1) + M
      return(make_out(scores, ballots = NULL, winner_idx = winners, picked = picked))
    }
    
    # No round-1 majority: choose top two finalists
    top_count <- max(counts1)
    top_set   <- which(counts1 == top_count)
    if (length(top_set) > 2) {
      finalists <- sample(top_set, 2)
    } else if (length(top_set) == 2) {
      finalists <- top_set
    } else {
      second_count <- sort(counts1, decreasing = TRUE)[2]
      runners <- which(counts1 == second_count)
      finalists <- c(top_set, sample(runners, 1))
    }
    a <- finalists[1]; b <- finalists[2]
    
    # Round 2: pairwise majority among finalists
    ab <- sum(utilities[, a] > utilities[, b], na.rm = TRUE)
    ba <- sum(utilities[, b] > utilities[, a], na.rm = TRUE)
    picked <- if (ab > ba) a else if (ba > ab) b else sample(c(a, b), 1)
    
    # Lexicographic scores: finalists above non-finalists; winner above runner-up; non-finalists by prop1
    M <- 2.0  # > any prop1 value to enforce lexicographic ordering
    scores <- prop1
    scores[finalists] <- max(prop1) + M            # both finalists above others
    scores[picked]    <- max(prop1) + M + 1.0      # winner above runner-up
    return(make_out(scores, ballots = NULL, winner_idx = finalists, picked = picked))
    
  } else if (scheme == "irv_r") {
    # Build strict order per voter (random tie-breaks)
    ranks <- t(apply(utilities, 1, function(x) order(-x, runif(length(x)))))
    
    active <- rep(TRUE, n_alternatives)
    elim_order <- integer(0)  # from first eliminated to winner
    
    repeat {
      # first active preference for each voter
      first_prefs <- apply(ranks, 1, function(ord) { for (idx in ord) if (active[idx]) return(idx); NA_integer_ })
      first_prefs <- first_prefs[!is.na(first_prefs)]
      counts <- tabulate(first_prefs, nbins = n_alternatives)
      
      total_active <- sum(counts)
      if (total_active == 0) {  # degenerate: pick random among active
        rem <- which(active); elim_order <- c(elim_order, sample(rem, length(rem))); break
      }
      
      # majority?
      if (max(counts) > total_active / 2) {
        winner <- which(counts == max(counts))
        winner <- sample(winner, 1)
        # append remaining active losers (except winner) in any order, then winner last
        rem <- which(active); rem <- setdiff(rem, winner)
        elim_order <- c(elim_order, rem, winner)
        break
      }
      
      # eliminate one with minimum first-preference count (break ties randomly)
      min_count <- min(counts[active])
      losers <- which(active & counts == min_count)
      drop <- sample(losers, 1)
      active[drop] <- FALSE
      elim_order <- c(elim_order, drop)
      
      if (sum(active) == 1) { elim_order <- c(elim_order, which(active)); break }
    }
    
    # Ranking score = reverse elimination order (later = larger)
    scores <- numeric(n_alternatives)
    scores[elim_order] <- seq_along(elim_order)   # winner gets max
    picked <- which.max(scores)
    return(make_out(scores, ballots = NULL, winner_idx = picked, picked = picked))
    
  } else if (scheme == "condorcet_r") {
    # Copeland completion (wins - losses). CW (if exists) will be unique max.
    wins <- matrix(0L, n_alternatives, n_alternatives)
    for (i in 1:(n_alternatives - 1)) for (j in (i + 1):n_alternatives) {
      ij <- sum(utilities[, i] > utilities[, j], na.rm = TRUE)
      ji <- sum(utilities[, j] > utilities[, i], na.rm = TRUE)
      if (ij > ji) wins[i, j] <- 1L else if (ji > ij) wins[j, i] <- 1L
    }
    w <- rowSums(wins); l <- colSums(wins)
    scores <- w - l
    
    # Condorcet winner check (beats all others)
    cw <- which(w == (n_alternatives - 1L))
    if (length(cw) == 1) {
      picked <- cw
      wset <- cw
    } else {
      best <- which(scores == max(scores))
      picked <- sample(best, 1)
      wset <- best
    }
    return(make_out(scores, ballots = NULL, winner_idx = wset, picked = picked))
    
  } else if (scheme == "typical_judgment_r") {
    # median + epsilon*(Pb - Pw)
    meds <- apply(utilities, 2, function(col) median(col, na.rm = TRUE))
    eps  <- 1e-3
    pb <- pw <- numeric(n_alternatives)
    for (j in 1:n_alternatives) {
      x <- utilities[, j]; x <- x[!is.na(x)]
      if (!length(x) || is.na(meds[j])) { pb[j] <- 0; pw[j] <- 0 } else {
        m <- meds[j]; n <- length(x)
        pb[j] <- sum(x > m) / n
        pw[j] <- sum(x < m) / n
      }
    }
    scores <- meds + eps * (pb - pw)
    best <- which(scores == max(scores, na.rm = TRUE))
    picked <- sample(best, 1)
    return(make_out(scores, ballots = NULL, winner_idx = best, picked = picked))
    
  } else if (scheme == "majority_judgment_r") {
    # median + epsilon*(Pb if Pb>Pw else -Pw)
    meds <- apply(utilities, 2, function(col) median(col, na.rm = TRUE))
    eps  <- 1e-3
    s2 <- numeric(n_alternatives)
    for (j in 1:n_alternatives) {
      x <- utilities[, j]; x <- x[!is.na(x)]
      if (!length(x) || is.na(meds[j])) { s2[j] <- 0 } else {
        m <- meds[j]; n <- length(x)
        pb <- sum(x > m) / n
        pw <- sum(x < m) / n
        s2[j] <- if (pb > pw) pb else if (pw > pb) -pw else 0
      }
    }
    scores <- meds + eps * s2
    best <- which(scores == max(scores, na.rm = TRUE))
    picked <- sample(best, 1)
    return(make_out(scores, ballots = NULL, winner_idx = best, picked = picked))
    } 
  else stop("Unsupported scheme: ", scheme)
}














