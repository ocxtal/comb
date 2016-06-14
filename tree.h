
/**
 * @file tree.h
 *
 * @brief an wrapper of ngx_rbtree.c in nginx (https://nginx.org/) core library
 */
#ifndef _TREE_H_INCLUDED
#define _TREE_H_INCLUDED

#include <stdint.h>


/**
 * @type rbtree_t
 */
typedef struct rbtree_s rbtree_t;

/**
 * @struct rbtree_node_s
 * @brief object must have a rbtree_node_t field at the head.
 */
struct rbtree_node_s {
	uint8_t pad[24];
	int64_t zero;				/* must be zeroed if external memory is used */
	int64_t key;
};
typedef struct rbtree_node_s rbtree_node_t;

/**
 * @struct rbtree_params_s
 */
struct rbtree_params_s {
	void *lmm;
};
typedef struct rbtree_params_s rbtree_params_t;
#define RBTREE_PARAMS(...)		( &((struct rbtree_params_s const) { __VA_ARGS__ }) )

/**
 * @fn rbtree_init
 */
rbtree_t *rbtree_init(uint64_t object_size, rbtree_params_t const *params);

/**
 * @fn rbtree_clean
 */
void rbtree_clean(rbtree_t *tree);

/**
 * @fn rbtree_flush
 */
void rbtree_flush(rbtree_t *tree);

/**
 * @fn rbtree_create_node
 * @brief create a new node (not inserted in the tree)
 */
rbtree_node_t *rbtree_create_node(rbtree_t *tree);

/**
 * @fn rbtree_insert
 * @brief insert a node
 */
void rbtree_insert(rbtree_t *tree, rbtree_node_t *node);

/**
 * @fn rbtree_remove
 * @brief remove a node, automatically freed if malloc'd with rbtree_reserve_node
 */
void rbtree_remove(rbtree_t *tree, rbtree_node_t *node);

/**
 * @fn rbtree_search_key
 * @brief search a node by key, returning the leftmost node
 */
rbtree_node_t *rbtree_search_key(rbtree_t *tree, int64_t key);

/**
 * @fn rbtree_search_key_left
 * @brief search a node by key. returns the nearest node in the left half of the tree if key was not found.
 */
rbtree_node_t *rbtree_search_key_left(rbtree_t *tree, int64_t key);

/**
 * @fn rbtree_search_key_right
 * @brief search a node by key. returns the nearest node in the right half of the tree if key was not found.
 */
rbtree_node_t *rbtree_search_key_right(rbtree_t *tree, int64_t key);

/**
 * @fn rbtree_left
 * @brief returns the left next node
 */
rbtree_node_t *rbtree_left(rbtree_t *tree, rbtree_node_t const *node);

/**
 * @fn rbtree_right
 * @brief returns the right next node
 */
rbtree_node_t *rbtree_right(rbtree_t *tree, rbtree_node_t const *node);

/**
 * @fn rbtree_walk
 * @breif iterate over tree
 */
typedef void (*rbtree_walk_t)(rbtree_node_t *node, void *ctx);
void rbtree_walk(rbtree_t *tree, rbtree_walk_t fn, void *ctx);




/* interval tree implementation */


/**
 * @type ivtree_t
 */
typedef struct rbtree_s ivtree_t;

/**
 * @struct ivtree_node_s
 */
struct ivtree_node_s {
	uint8_t pad[24];
	int64_t zero;				/* must be zeroed if external memory is used */
	int64_t lkey;
	int64_t rkey;
	int64_t reserved;
};
typedef struct ivtree_node_s ivtree_node_t;

/**
 * @type ivtree_iter_t
 */
typedef struct ivtree_iter_s ivtree_iter_t;

/**
 * @type ivtree_params_s
 */
typedef struct rvtree_params_s ivtree_params_t;
#define IVTREE_PARAMS(...)		( &((struct rbtree_params_s const) { __VA_ARGS__ }) )

/**
 * @fn ivtree_init
 */
ivtree_t *ivtree_init(uint64_t object_size, ivtree_params_t const *params);

/**
 * @fn ivtree_clean
 */
void ivtree_clean(ivtree_t *tree);

/**
 * @fn ivtree_flush
 */
void ivtree_flush(ivtree_t *tree);

/**
 * @fn ivtree_create_node
 * @brief create a new node (not inserted in the tree)
 */
ivtree_node_t *ivtree_create_node(ivtree_t *tree);

/**
 * @fn ivtree_insert
 * @brief insert a node
 */
void ivtree_insert(ivtree_t *tree, ivtree_node_t *node);

/**
 * @fn ivtree_remove
 * @brief remove a node, automatically freed if malloc'd with ivtree_reserve_node
 */
void ivtree_remove(ivtree_t *tree, ivtree_node_t *node);

/**
 * @fn ivtree_contained
 * @brief return a set of sections contained in [lkey, rkey)
 */
ivtree_iter_t *ivtree_contained(ivtree_t *tree, int64_t lkey, int64_t rkey);

/**
 * @fn ivtree_containing
 * @brief return a set of sections containing [lkey, rkey)
 */
ivtree_iter_t *ivtree_containing(ivtree_t *tree, int64_t lkey, int64_t rkey);

/**
 * @fn ivtree_intersect
 * @brief return a set of sections intersect with [lkey, rkey)
 */
ivtree_iter_t *ivtree_intersect(ivtree_t *tree, int64_t lkey, int64_t rkey);

/**
 * @fn ivtree_next
 */
ivtree_node_t *ivtree_next(ivtree_iter_t *iter);

/**
 * @fn ivtree_iter_clean
 */
void ivtree_iter_clean(ivtree_iter_t *iter);

/**
 * @fn ivtree_update_node
 */
void ivtree_update_node(ivtree_t *tree, ivtree_node_t *node);

/**
 * @fn ivtree_walk
 * @breif iterate over tree
 */
typedef void (*ivtree_walk_t)(ivtree_node_t *node, void *ctx);
void ivtree_walk(ivtree_t *tree, ivtree_walk_t fn, void *ctx);


#endif
/**
 * end of tree.h
 */
