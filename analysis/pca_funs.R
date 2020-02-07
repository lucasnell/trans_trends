library(tidyverse)

# predicted values from model
pred_fn <- function(d_, beta_, int_) {
    if ("dist" %in% beta_$coef) {
        dd <- sym("dist")
        if (!"dist" %in% colnames(d_)) {
            stop("colnames in d_ and beta_$coef don't agree")
        }
    } else if ("distance" %in% beta_$coef) {
        dd <- sym("distance")
        if (!"distance" %in% colnames(d_)) {
            stop("colnames in d_ and beta_$coef don't agree")
        }
    } else stop("Neither dist nor distance found in `pred_fn`")
    d_ %>%
        expand(taxon,
               midges_z = seq(min(midges_z), max(midges_z), 1),
               time_z = seq(min(time_z), max(time_z), 1),
               dist_z = seq(min(dist_z), max(dist_z), 1)) %>%
        full_join(beta_ %>%
                      select(taxon, coef, mi) %>%
                      spread(coef, mi)) %>%
        full_join(int_ %>% select(taxon, mi) %>% rename(int = mi)) %>%
        mutate(y = int + midges*midges_z + time*time_z + !!dd*dist_z) %>%
        group_by(taxon) %>%
        mutate(id = row_number()) %>%
        ungroup()

}

# pca on predicted values
pca_fn <- function(pred_) {
    prcomp(pred_ %>%
               select(id, taxon, y) %>%
               spread(taxon, y) %>%
               select(-id) %>%
               as.matrix(),
           center = F, scale= F)
}


# axes from pca
axes_fn <- function(pred_, pca_) {
    pred_ %>%
        select(id, midges_z, time_z, dist_z) %>%
        full_join(as_tibble(pca_$x) %>%
                      mutate(id = row_number()))
}

# taxon vectors from pca
taxon_vec_fn <- function(d_, pca_) {
    taxa_short_ <- tibble(taxon = d_$taxon %>% unique()) %>%
        mutate(id = row_number())

    pca_$rotation %>%
        as_tibble() %>%
        mutate(id = row_number()) %>%
        left_join(taxa_short_) %>%
        select(taxon, matches("PC"))
}

# rotation of observed data
obs_rot_fn <- function(d_, pca_) {
    obs_ <- d_ %>%
        group_by(taxon) %>%
        mutate(id = row_number()) %>%
        ungroup()

    obs_ %>%
        select(id, plot, trans, midges_z, time_z, dist_z) %>%
        full_join(as_tibble(predict(pca_, obs_ %>%
                                        select(id, taxon, y) %>%
                                        spread(taxon, y) %>%
                                        select(-id) %>%
                                        as.matrix())) %>%
                      mutate(id = row_number()))
}

# variance explained in observed data by rotation
obs_exp_fn <- function(obs_rot_) {
    obs_rot_ %>% select(matches("PC")) %>%
        apply(2, sd) %>%
        {.^2/sum(.^2)} %>%
        rbind(.,cumsum(.)) %>%
        as_tibble() %>%
        mutate(type = c("individual","cumulative")) %>%
        select(type, matches("PC"))
}

# run pca and package
pred_pca_fn <- function(d_, beta_, int_) {
    pred_ <- pred_fn(d_, beta_, int_)
    pca_ <- pca_fn(pred_)
    axes_ <- axes_fn(pred_, pca_)
    taxon_vec_ <- taxon_vec_fn(d_, pca_)
    obs_rot_ <- obs_rot_fn(d_, pca_)
    obs_exp_ <- obs_exp_fn(obs_rot_)

    return(list(pred = pred_, pca = pca_, axes = axes_, taxon_vec = taxon_vec_,
                obs_rot = obs_rot_, obs_exp = obs_exp_))
}
