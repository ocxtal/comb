
/*
 * Copyright (C) Igor Sysoev
 * Copyright (C) Nginx, Inc.
 */


#ifndef _NGX_RBTREE_H_INCLUDED_
#define _NGX_RBTREE_H_INCLUDED_


// #include <ngx_config.h>
// #include <ngx_core.h>

#include <stdint.h>         // for uint64_t

// typedef ngx_uint_t  ngx_rbtree_key_t;
// typedef ngx_int_t   ngx_rbtree_key_int_t;

/**
 * modified to hold 64bit key-value pairs
 */
// typedef int64_t ngx_rbtree_key_t;


typedef struct ngx_rbtree_node_s  ngx_rbtree_node_t;

struct ngx_rbtree_node_s {
    ngx_rbtree_node_t       *parent;
    ngx_rbtree_node_t       *left;
    ngx_rbtree_node_t       *right;
    uint8_t                 color;
    uint8_t                 data;
    uint8_t                 pad[6];
    int64_t                 key;
};


typedef struct ngx_rbtree_s  ngx_rbtree_t;

typedef void (*ngx_rbtree_insert_pt) (ngx_rbtree_node_t *root,
    ngx_rbtree_node_t *node, ngx_rbtree_node_t *sentinel);

struct ngx_rbtree_s {
    ngx_rbtree_node_t     *root;
    ngx_rbtree_node_t     *sentinel;
    // ngx_rbtree_insert_pt   insert;
};


#define ngx_rbtree_init(tree, s, i)                                           \
    ngx_rbtree_sentinel_init(s);                                              \
    (tree)->root = s;                                                         \
    (tree)->sentinel = s;


void ngx_rbtree_insert(ngx_rbtree_t *tree, ngx_rbtree_node_t *node);
void ngx_rbtree_delete(ngx_rbtree_t *tree, ngx_rbtree_node_t *node);

/*
 * search functions
 * find_key return the leftmost node
 * added 2015/11/06
 */
ngx_rbtree_node_t *ngx_rbtree_find_key(ngx_rbtree_t *tree, int64_t key);
ngx_rbtree_node_t *ngx_rbtree_find_key_left(ngx_rbtree_t *tree, int64_t key);
ngx_rbtree_node_t *ngx_rbtree_find_key_right(ngx_rbtree_t *tree, int64_t key);
ngx_rbtree_node_t *ngx_rbtree_find_left(ngx_rbtree_t *tree, ngx_rbtree_node_t *node);
ngx_rbtree_node_t *ngx_rbtree_find_right(ngx_rbtree_t *tree, ngx_rbtree_node_t *node);

typedef void (*ngx_rbtree_walk_pt) (ngx_rbtree_node_t **node, ngx_rbtree_node_t *sentinel, void *ctx);
void ngx_rbtree_walk(ngx_rbtree_t *tree, ngx_rbtree_walk_pt walk, void *ctx);

#define ngx_rbt_red(node)               ((node)->color = 1)
#define ngx_rbt_black(node)             ((node)->color = 0)
#define ngx_rbt_is_red(node)            ((node)->color)
#define ngx_rbt_is_black(node)          (!ngx_rbt_is_red(node))
#define ngx_rbt_copy_color(n1, n2)      (n1->color = n2->color)


/* a sentinel must be black */

#define ngx_rbtree_sentinel_init(node)  ngx_rbt_black(node)


/* interval tree (augumented tree) */

typedef struct ngx_ivtree_node_s ngx_ivtree_node_t;

struct ngx_ivtree_node_s {
    ngx_ivtree_node_t       *parent;
    ngx_ivtree_node_t       *left;
    ngx_ivtree_node_t       *right;
    uint8_t                 color;
    uint8_t                 data;
    uint8_t                 pad[6];
    int64_t                 lkey;
    int64_t                 rkey;
    int64_t                 rkey_max;
};


typedef struct ngx_ivtree_s  ngx_ivtree_t;


struct ngx_ivtree_s {
    ngx_ivtree_node_t     *root;
    ngx_ivtree_node_t     *sentinel;
    // ngx_rbtree_insert_pt   insert;
};


#define ngx_ivtree_init(tree, s, i)     ngx_rbtree_init(tree, s, i)


void ngx_ivtree_insert(ngx_ivtree_t *tree, ngx_ivtree_node_t *node);
void ngx_ivtree_delete(ngx_ivtree_t *tree, ngx_ivtree_node_t *node);


#endif /* _NGX_RBTREE_H_INCLUDED_ */
