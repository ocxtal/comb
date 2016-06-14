
/*
 * Copyright (C) Igor Sysoev
 * Copyright (C) Nginx, Inc.
 */


// #include <ngx_config.h>
// #include <ngx_core.h>
#include <stdlib.h>
#include "ngx_rbtree.h"
#include "log.h"

/*
 * The red-black tree code is based on the algorithm described in
 * the "Introduction to Algorithms" by Cormen, Leiserson and Rivest.
 */

static inline void ngx_rbtree_insert_value(ngx_rbtree_node_t *temp,
    ngx_rbtree_node_t *node, ngx_rbtree_node_t *sentinel);
static inline void ngx_rbtree_left_rotate(ngx_rbtree_node_t **root,
    ngx_rbtree_node_t *sentinel, ngx_rbtree_node_t *node);
static inline void ngx_rbtree_right_rotate(ngx_rbtree_node_t **root,
    ngx_rbtree_node_t *sentinel, ngx_rbtree_node_t *node);


static inline ngx_rbtree_node_t *
ngx_rbtree_min(ngx_rbtree_node_t *node, ngx_rbtree_node_t *sentinel)
{
    debug("node(%p)", node);
    while (node->left != sentinel) {
        debug("node(%p)", node);
        node = node->left;
    }

    return node;
}


void
ngx_rbtree_insert(ngx_rbtree_t *tree, ngx_rbtree_node_t *node)
{
    ngx_rbtree_node_t  **root, *temp, *sentinel;

    /* a binary tree insert */

    root = (ngx_rbtree_node_t **) &tree->root;
    sentinel = tree->sentinel;

    if (*root == sentinel) {
        node->parent = NULL;
        node->left = sentinel;
        node->right = sentinel;
        ngx_rbt_black(node);
        *root = node;

        return;
    }

    ngx_rbtree_insert_value(*root, node, sentinel);

    /* re-balance tree */

    while (node != *root && ngx_rbt_is_red(node->parent)) {

        if (node->parent == node->parent->parent->left) {
            temp = node->parent->parent->right;

            if (ngx_rbt_is_red(temp)) {
                ngx_rbt_black(node->parent);
                ngx_rbt_black(temp);
                ngx_rbt_red(node->parent->parent);
                node = node->parent->parent;

            } else {
                if (node == node->parent->right) {
                    node = node->parent;
                    ngx_rbtree_left_rotate(root, sentinel, node);
                }

                ngx_rbt_black(node->parent);
                ngx_rbt_red(node->parent->parent);
                ngx_rbtree_right_rotate(root, sentinel, node->parent->parent);
            }

        } else {
            temp = node->parent->parent->left;

            if (ngx_rbt_is_red(temp)) {
                ngx_rbt_black(node->parent);
                ngx_rbt_black(temp);
                ngx_rbt_red(node->parent->parent);
                node = node->parent->parent;

            } else {
                if (node == node->parent->left) {
                    node = node->parent;
                    ngx_rbtree_right_rotate(root, sentinel, node);
                }

                ngx_rbt_black(node->parent);
                ngx_rbt_red(node->parent->parent);
                ngx_rbtree_left_rotate(root, sentinel, node->parent->parent);
            }
        }
    }

    ngx_rbt_black(*root);
}


static inline void
ngx_rbtree_insert_value(ngx_rbtree_node_t *temp, ngx_rbtree_node_t *node,
    ngx_rbtree_node_t *sentinel)
{
    ngx_rbtree_node_t  **p;

    for ( ;; ) {

        p = (node->key < temp->key) ? &temp->left : &temp->right;

        if (*p == sentinel) {
            break;
        }

        temp = *p;
    }

    *p = node;
    node->parent = temp;
    node->left = sentinel;
    node->right = sentinel;
    ngx_rbt_red(node);
}


void
ngx_rbtree_delete(ngx_rbtree_t *tree, ngx_rbtree_node_t *node)
{
    uint8_t           red;
    ngx_rbtree_node_t  **root, *sentinel, *subst, *temp, *w;

    /* a binary tree delete */

    root = (ngx_rbtree_node_t **) &tree->root;
    sentinel = tree->sentinel;

    if (node->left == sentinel) {
        temp = node->right;
        subst = node;

    } else if (node->right == sentinel) {
        temp = node->left;
        subst = node;

    } else {
        subst = ngx_rbtree_min(node->right, sentinel);
        if (subst->left != sentinel) {
            temp = subst->left;
        } else {
            temp = subst->right;
        }

        // temp = subst->right;
    }

    if (subst == *root) {
        *root = temp;
        ngx_rbt_black(temp);

        /* DEBUG stuff */
        node->left = NULL;
        node->right = NULL;
        node->parent = NULL;
        node->key = 0;

        return;
    }

    red = ngx_rbt_is_red(subst);

    if (subst == subst->parent->left) {
        subst->parent->left = temp;

    } else {
        subst->parent->right = temp;
    }

    if (subst == node) {

        temp->parent = subst->parent;

    } else {

        if (subst->parent == node) {
            temp->parent = subst;

        } else {
            temp->parent = subst->parent;
        }

        subst->left = node->left;
        subst->right = node->right;
        subst->parent = node->parent;
        ngx_rbt_copy_color(subst, node);

        if (node == *root) {
            *root = subst;

        } else {
            if (node == node->parent->left) {
                node->parent->left = subst;
            } else {
                node->parent->right = subst;
            }
        }

        if (subst->left != sentinel) {
            subst->left->parent = subst;
        }

        if (subst->right != sentinel) {
            subst->right->parent = subst;
        }
    }

    /* DEBUG stuff */
    node->left = NULL;
    node->right = NULL;
    node->parent = NULL;
    node->key = 0;

    if (red) {
        return;
    }

    /* a delete fixup */

    while (temp != *root && ngx_rbt_is_black(temp)) {

        if (temp == temp->parent->left) {
            w = temp->parent->right;

            if (ngx_rbt_is_red(w)) {
                ngx_rbt_black(w);
                ngx_rbt_red(temp->parent);
                ngx_rbtree_left_rotate(root, sentinel, temp->parent);
                w = temp->parent->right;
            }

            if (ngx_rbt_is_black(w->left) && ngx_rbt_is_black(w->right)) {
                ngx_rbt_red(w);
                temp = temp->parent;

            } else {
                if (ngx_rbt_is_black(w->right)) {
                    ngx_rbt_black(w->left);
                    ngx_rbt_red(w);
                    ngx_rbtree_right_rotate(root, sentinel, w);
                    w = temp->parent->right;
                }

                ngx_rbt_copy_color(w, temp->parent);
                ngx_rbt_black(temp->parent);
                ngx_rbt_black(w->right);
                ngx_rbtree_left_rotate(root, sentinel, temp->parent);
                temp = *root;
            }

        } else {
            w = temp->parent->left;

            if (ngx_rbt_is_red(w)) {
                ngx_rbt_black(w);
                ngx_rbt_red(temp->parent);
                ngx_rbtree_right_rotate(root, sentinel, temp->parent);
                w = temp->parent->left;
            }

            if (ngx_rbt_is_black(w->left) && ngx_rbt_is_black(w->right)) {
                ngx_rbt_red(w);
                temp = temp->parent;

            } else {
                if (ngx_rbt_is_black(w->left)) {
                    ngx_rbt_black(w->right);
                    ngx_rbt_red(w);
                    ngx_rbtree_left_rotate(root, sentinel, w);
                    w = temp->parent->left;
                }

                ngx_rbt_copy_color(w, temp->parent);
                ngx_rbt_black(temp->parent);
                ngx_rbt_black(w->left);
                ngx_rbtree_right_rotate(root, sentinel, temp->parent);
                temp = *root;
            }
        }
    }

    ngx_rbt_black(temp);
}


ngx_rbtree_node_t *
ngx_rbtree_find_key(ngx_rbtree_t *tree, int64_t key)
{
    ngx_rbtree_node_t *node = tree->root;
    ngx_rbtree_node_t *sentinel = tree->sentinel;

    while(node != sentinel) {
        if(key < node->key) {
            node = node->left;
        } else if(key > node->key) {
            node = node->right;
        } else {
            /* key == node->key */
            while(node->left != sentinel && key == node->left->key) {
                node = node->left;
            }
            return(node);
        }
    }
    return(NULL);
}


ngx_rbtree_node_t *
ngx_rbtree_find_key_right(ngx_rbtree_t *tree, int64_t key)
{
    ngx_rbtree_node_t *node = tree->root;
    ngx_rbtree_node_t *sentinel = tree->sentinel;

    if(node == sentinel) {
        return(NULL);
    }
    for(;;) {
        if(key < node->key) {
            if(node->left == sentinel) {
                return(node);
            }
            node = node->left;
        } else if(key > node->key) {
            if(node->right == sentinel) {
                return(ngx_rbtree_find_right(tree, node));
            }
            node = node->right;
        } else {
            /* key == node->key */
            while(node->left != sentinel && key == node->left->key) {
                node = node->left;
            }
            return(node);
        }
    }
    return(NULL);
}


ngx_rbtree_node_t *
ngx_rbtree_find_key_left(ngx_rbtree_t *tree, int64_t key)
{
    ngx_rbtree_node_t *node = tree->root;
    ngx_rbtree_node_t *sentinel = tree->sentinel;

    if(node == sentinel) {
        return(NULL);
    }
    for(;;) {
        if(key < node->key) {
            if(node->left == sentinel) {
                return(ngx_rbtree_find_left(tree, node));
            }
            node = node->left;
        } else if(key > node->key) {
            if(node->right == sentinel) {
                return(node);
            }
            node = node->right;
        } else {
            /* key == node->key */
            while(node->left != sentinel && key == node->left->key) {
                node = node->left;
            }
            return(node);
        }
    }
    return(NULL);
}


ngx_rbtree_node_t *
ngx_rbtree_find_right(ngx_rbtree_t *tree, ngx_rbtree_node_t *node)
{
    ngx_rbtree_node_t *sentinel = tree->sentinel;

    if(node == sentinel) {
        return(NULL);
    }

    if(node->right != sentinel) {
        node = node->right;
        while(node->left != sentinel) {
            node = node->left;
        }
        return(node);
    }

    /* when the node has no right child */
    while(node->parent != NULL && node == node->parent->right) {
        node = node->parent;
    }
    return(node->parent);
}


ngx_rbtree_node_t *
ngx_rbtree_find_left(ngx_rbtree_t *tree, ngx_rbtree_node_t *node)
{
    ngx_rbtree_node_t *sentinel = tree->sentinel;

    if(node == sentinel) {
        return(NULL);
    }

    if(node->left != sentinel) {
        node = node->left;
        while(node->right != sentinel) {
            node = node->right;
        }
        return(node);
    }

    /* when the node has no left child */
    while(node->parent != NULL && node == node->parent->left) {
        node = node->parent;
    }
    return(node->parent);
}

static void
ngx_rbtree_walk_intl(ngx_rbtree_node_t *node, ngx_rbtree_node_t *sentinel, ngx_rbtree_walk_pt walk, void *ctx)
{

    if(node->left != sentinel) {
        ngx_rbtree_walk_intl(node->left, sentinel, walk, ctx);
    }
    if(node->right != sentinel) {
        ngx_rbtree_walk_intl(node->right, sentinel, walk, ctx);
    }
    walk(&node, sentinel, ctx);
}

void
ngx_rbtree_walk(ngx_rbtree_t *tree, ngx_rbtree_walk_pt walk, void *ctx)
{
    ngx_rbtree_node_t *node = tree->root;
    ngx_rbtree_node_t *sentinel = tree->sentinel;

    if(node != sentinel) {
        ngx_rbtree_walk_intl(node, sentinel, walk, ctx);
    }

    return;
}


static inline void
ngx_rbtree_left_rotate(ngx_rbtree_node_t **root, ngx_rbtree_node_t *sentinel,
    ngx_rbtree_node_t *node)
{
    ngx_rbtree_node_t  *temp;

    temp = node->right;
    node->right = temp->left;

    if (temp->left != sentinel) {
        temp->left->parent = node;
    }

    temp->parent = node->parent;

    if (node == *root) {
        *root = temp;

    } else if (node == node->parent->left) {
        node->parent->left = temp;

    } else {
        node->parent->right = temp;
    }

    temp->left = node;
    node->parent = temp;
}


static inline void
ngx_rbtree_right_rotate(ngx_rbtree_node_t **root, ngx_rbtree_node_t *sentinel,
    ngx_rbtree_node_t *node)
{
    ngx_rbtree_node_t  *temp;

    temp = node->left;
    node->left = temp->right;

    if (temp->right != sentinel) {
        temp->right->parent = node;
    }

    temp->parent = node->parent;

    if (node == *root) {
        *root = temp;

    } else if (node == node->parent->right) {
        node->parent->right = temp;

    } else {
        node->parent->left = temp;
    }

    temp->right = node;
    node->parent = temp;
}

#if 1
/* interval tree implementation */
/* max and min */
#define MAX2(x,y)       ( (x) > (y) ? (x) : (y) )
#define MAX3(x,y,z)     ( MAX2(x, MAX2(y, z)) )

#define MIN2(x,y)       ( (x) < (y) ? (x) : (y) )
#define MIN3(x,y,z)     ( MIN2(x, MIN2(y, z)) )


static inline void ngx_ivtree_left_rotate(ngx_ivtree_node_t **root,
    ngx_ivtree_node_t *sentinel, ngx_ivtree_node_t *node);
static inline void ngx_ivtree_right_rotate(ngx_ivtree_node_t **root,
    ngx_ivtree_node_t *sentinel, ngx_ivtree_node_t *node);

void
ngx_ivtree_update_key(ngx_ivtree_t *tree, ngx_ivtree_node_t *node)
{
    while(node->parent != NULL) {
        debug("parent(%p, %lld, %lld, %lld), node(%p, %lld, %lld, %lld)",
            node->parent, node->parent->lkey, node->parent->rkey, node->parent->rkey_max,
            node, node->lkey, node->rkey, node->rkey_max);
        if(node->parent->rkey_max < node->rkey_max) {
            node->parent->rkey_max = node->rkey_max;
            break;
        }
        node = node->parent;
    }
    return;
}

void
ngx_ivtree_insert(ngx_ivtree_t *tree, ngx_ivtree_node_t *node)
{
    ngx_ivtree_node_t  **root, *temp, *sentinel;

    /* a binary tree insert */

    root = (ngx_ivtree_node_t **) &tree->root;
    sentinel = tree->sentinel;

    if (*root == sentinel) {
        node->parent = NULL;
        node->left = sentinel;
        node->right = sentinel;
        node->rkey_max = node->rkey;
        ngx_rbt_black(node);
        *root = node;

        return;
    }

    debug("insert value node(%p, %lld, %lld)", node, node->lkey, node->rkey);
    ngx_rbtree_insert_value(
        (ngx_rbtree_node_t *)*root,
        (ngx_rbtree_node_t *)node,
        (ngx_rbtree_node_t *)sentinel);
    node->rkey_max = node->rkey;      /* inital max */
    ngx_ivtree_update_key(tree, node);

    /* re-balance tree */

    while (node != *root && ngx_rbt_is_red(node->parent)) {

        if (node->parent == node->parent->parent->left) {
            temp = node->parent->parent->right;

            if (ngx_rbt_is_red(temp)) {
                ngx_rbt_black(node->parent);
                ngx_rbt_black(temp);
                ngx_rbt_red(node->parent->parent);
                node = node->parent->parent;

            } else {
                if (node == node->parent->right) {
                    node = node->parent;
                    ngx_ivtree_left_rotate(root, sentinel, node);
                }

                ngx_rbt_black(node->parent);
                ngx_rbt_red(node->parent->parent);
                ngx_ivtree_right_rotate(root, sentinel, node->parent->parent);
            }

        } else {
            temp = node->parent->parent->left;

            if (ngx_rbt_is_red(temp)) {
                ngx_rbt_black(node->parent);
                ngx_rbt_black(temp);
                ngx_rbt_red(node->parent->parent);
                node = node->parent->parent;

            } else {
                if (node == node->parent->left) {
                    node = node->parent;
                    ngx_ivtree_right_rotate(root, sentinel, node);
                }

                ngx_rbt_black(node->parent);
                ngx_rbt_red(node->parent->parent);
                ngx_ivtree_left_rotate(root, sentinel, node->parent->parent);
            }
        }
    }

    ngx_rbt_black(*root);
}


void
ngx_ivtree_delete(ngx_ivtree_t *tree, ngx_ivtree_node_t *node)
{
    uint8_t           red;
    ngx_ivtree_node_t  **root, *sentinel, *subst, *temp, *w;

    /* a binary tree delete */

    root = (ngx_ivtree_node_t **) &tree->root;
    sentinel = tree->sentinel;

    if (node->left == sentinel) {
        temp = node->right;
        subst = node;

    } else if (node->right == sentinel) {
        temp = node->left;
        subst = node;

    } else {
        subst = (ngx_ivtree_node_t *)ngx_rbtree_min(
            (ngx_rbtree_node_t *)node->right,
            (ngx_rbtree_node_t *)sentinel);
        if (subst->left != sentinel) {
            temp = subst->left;
        } else {
            temp = subst->right;
        }

        // temp = subst->right;
    }

    if (subst == *root) {
        *root = temp;
        ngx_rbt_black(temp);

        /* DEBUG stuff */
        node->left = NULL;
        node->right = NULL;
        node->parent = NULL;
        node->lkey = 0;
        node->rkey = 0;

        return;
    }

    red = ngx_rbt_is_red(subst);

    if (subst == subst->parent->left) {
        subst->parent->left = temp;

    } else {
        subst->parent->right = temp;
    }
    subst->parent->rkey_max = MAX3(
        subst->parent->rkey,
        subst->parent->left->rkey_max,
        subst->parent->right->rkey_max);

    if (subst == node) {

        temp->parent = subst->parent;

    } else {

        if (subst->parent == node) {
            temp->parent = subst;

        } else {
            temp->parent = subst->parent;
        }

        subst->left = node->left;
        subst->right = node->right;
        subst->parent = node->parent;
        ngx_rbt_copy_color(subst, node);
        subst->rkey_max = MAX3(
            subst->rkey,
            subst->left->rkey_max,
            subst->right->rkey_max);

        if (node == *root) {
            *root = subst;

        } else {
            if (node == node->parent->left) {
                node->parent->left = subst;
            } else {
                node->parent->right = subst;
            }
            node->parent->rkey_max = MAX3(
                node->parent->rkey,
                node->parent->left->rkey_max,
                node->parent->right->rkey_max);

        }

        if (subst->left != sentinel) {
            subst->left->parent = subst;
        }

        if (subst->right != sentinel) {
            subst->right->parent = subst;
        }
    }

    /* DEBUG stuff */
    node->left = NULL;
    node->right = NULL;
    node->parent = NULL;
    node->lkey = 0;
    node->rkey = 0;

    if (red) {
        return;
    }

    /* a delete fixup */
    while (temp != *root && ngx_rbt_is_black(temp)) {

        if (temp == temp->parent->left) {
            w = temp->parent->right;

            if (ngx_rbt_is_red(w)) {
                ngx_rbt_black(w);
                ngx_rbt_red(temp->parent);
                ngx_ivtree_left_rotate(root, sentinel, temp->parent);
                w = temp->parent->right;
            }

            if (ngx_rbt_is_black(w->left) && ngx_rbt_is_black(w->right)) {
                ngx_rbt_red(w);
                temp = temp->parent;

            } else {
                if (ngx_rbt_is_black(w->right)) {
                    ngx_rbt_black(w->left);
                    ngx_rbt_red(w);
                    ngx_ivtree_right_rotate(root, sentinel, w);
                    w = temp->parent->right;
                }

                ngx_rbt_copy_color(w, temp->parent);
                ngx_rbt_black(temp->parent);
                ngx_rbt_black(w->right);
                ngx_ivtree_left_rotate(root, sentinel, temp->parent);
                temp = *root;
            }

        } else {
            w = temp->parent->left;

            if (ngx_rbt_is_red(w)) {
                ngx_rbt_black(w);
                ngx_rbt_red(temp->parent);
                ngx_ivtree_right_rotate(root, sentinel, temp->parent);
                w = temp->parent->left;
            }

            if (ngx_rbt_is_black(w->left) && ngx_rbt_is_black(w->right)) {
                ngx_rbt_red(w);
                temp = temp->parent;

            } else {
                if (ngx_rbt_is_black(w->left)) {
                    ngx_rbt_black(w->right);
                    ngx_rbt_red(w);
                    ngx_ivtree_left_rotate(root, sentinel, w);
                    w = temp->parent->left;
                }

                ngx_rbt_copy_color(w, temp->parent);
                ngx_rbt_black(temp->parent);
                ngx_rbt_black(w->left);
                ngx_ivtree_right_rotate(root, sentinel, temp->parent);
                temp = *root;
            }
        }
    }

    ngx_rbt_black(temp);
}


static inline void
ngx_ivtree_left_rotate(ngx_ivtree_node_t **root, ngx_ivtree_node_t *sentinel,
    ngx_ivtree_node_t *node)
{
    ngx_ivtree_node_t  *temp;

    temp = node->right;
    node->right = temp->left;

    if (temp->left != sentinel) {
        temp->left->parent = node;
    }

    temp->parent = node->parent;

    if (node == *root) {
        *root = temp;

    } else if (node == node->parent->left) {
        node->parent->left = temp;

    } else {
        node->parent->right = temp;
    }

    temp->left = node;
    node->parent = temp;

    node->rkey_max = MAX3(
        node->rkey,
        node->left->rkey_max,
        node->right->rkey_max);
    node->parent->rkey_max = MAX3(
        node->parent->rkey,
        node->parent->left->rkey_max,
        node->parent->right->rkey_max);

}


static inline void
ngx_ivtree_right_rotate(ngx_ivtree_node_t **root, ngx_ivtree_node_t *sentinel,
    ngx_ivtree_node_t *node)
{
    ngx_ivtree_node_t  *temp;

    temp = node->left;
    node->left = temp->right;

    if (temp->right != sentinel) {
        temp->right->parent = node;
    }

    temp->parent = node->parent;

    if (node == *root) {
        *root = temp;

    } else if (node == node->parent->right) {
        node->parent->right = temp;

    } else {
        node->parent->left = temp;
    }

    temp->right = node;
    node->parent = temp;

    node->rkey_max = MAX3(
        node->rkey,
        node->left->rkey_max,
        node->right->rkey_max);
    node->parent->rkey_max = MAX3(
        node->parent->rkey,
        node->parent->left->rkey_max,
        node->parent->right->rkey_max);

}

#endif