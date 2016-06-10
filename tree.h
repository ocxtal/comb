
/**
 * @file tree.h
 *
 * @brief an wrapper of ngx_rbtree.c in nginx (https://nginx.org/) core library
 */
#ifndef _TREE_H_INCLUDED
#define _TREE_H_INCLUDED

#include <stdint.h>


/**
 * @type tree_t
 */
typedef struct tree_s tree_t;

/**
 * @struct tree_node_s
 */
struct tree_node_s {
	uint8_t pad[24];
	int64_t zero;				/* must be zeroed if external memory is used */
	int64_t key;
};
typedef struct tree_node_s tree_node_t;

/**
 * @struct tree_params_s
 */
struct tree_params_s {
	void *lmm;
};
typedef struct tree_params_s tree_params_t;
#define TREE_PARAMS(...)		( &((struct tree_params_s const) { __VA_ARGS__ }) )

/**
 * @macro tree_get_object
 * @brief retrieve object pointer from node pointer
 */
#define TREE_OBJECT_OFFSET		( 32 )
#define tree_get_object(x)		( (void *)((uint8_t *)(x) + TREE_OBJECT_OFFSET) )


/**
 * @fn tree_init
 */
tree_t *tree_init(uint64_t object_size, tree_params_t const *params);

/**
 * @fn tree_clean
 */
void tree_clean(tree_t *tree);

/**
 * @fn tree_flush
 */
void tree_flush(tree_t *tree);

/**
 * @fn tree_create_node
 * @brief create a new node (not inserted in the tree)
 */
tree_node_t *tree_create_node(tree_t *tree);

/**
 * @fn tree_insert
 * @brief insert a node
 */
void tree_insert(tree_t *tree, tree_node_t *node);

/**
 * @fn tree_remove
 * @brief remove a node, automatically freed if malloc'd with tree_reserve_node
 */
void tree_remove(tree_t *tree, tree_node_t *node);

/**
 * @fn tree_search_key
 * @brief search a node by key, returning the leftmost node
 */
tree_node_t *tree_search_key(tree_t *tree, int64_t key);

/**
 * @fn tree_search_key_left
 * @brief search a node by key. returns the nearest node in the left half of the tree if key was not found.
 */
tree_node_t *tree_search_key_left(tree_t *tree, int64_t key);

/**
 * @fn tree_search_key_right
 * @brief search a node by key. returns the nearest node in the right half of the tree if key was not found.
 */
tree_node_t *tree_search_key_right(tree_t *tree, int64_t key);

/**
 * @fn tree_left
 * @brief returns the left next node
 */
tree_node_t *tree_left(tree_t *tree, tree_node_t const *node);

/**
 * @fn tree_right
 * @brief returns the right next node
 */
tree_node_t *tree_right(tree_t *tree, tree_node_t const *node);

/**
 * @fn tree_walk
 * @breif iterate over tree
 */
typedef void (*tree_walk_t)(tree_node_t *node, void *ctx);
void tree_walk(tree_t *tree, tree_walk_t fn, void *ctx);

/**
 * end of tree.c
 */


#endif
/**
 * end of tree.h
 */
