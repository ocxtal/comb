
/**
 * @file tree.c
 *
 * @brief an wrapper of ngx_rbtree.c in nginx (https://nginx.org/) core library
 */

#define UNITTEST_UNIQUE_ID		59
#include "unittest.h"

#include <stdint.h>
#include <stdlib.h>
#include "ngx_rbtree.h"
#include "lmm.h"
#include "log.h"
#include "sassert.h"
#include "tree.h"


/* constants */
#define  TREE_INIT_ELEM_CNT			( 64 )

/* roundup */
#define _roundup(x, base)			( ((x) + (base) - 1) & ~((base) - 1) )

/**
 * @struct tree_vec_s
 */
struct tree_vec_s {
	struct tree_vec_s *next;
	struct tree_vec_s *prev;
	int64_t cnt;
	int64_t pad;
	uint8_t v[];
};

/**
 * @struct tree_s
 */
struct tree_s {
	lmm_t *lmm;
	uint32_t object_size;
	uint32_t pad;
	struct tree_params_s params;

	/* vector pointers */
	uint8_t *v;
	int64_t vrem;
	struct tree_vec_s *vhead, *vroot;

	/* tree */
	ngx_rbtree_t t;
	ngx_rbtree_node_t sentinel;

	/* working buffer */
	tree_walk_t fn;
	void *ctx;
};

/**
 * @struct tree_node_intl_s
 */
struct tree_node_intl_s {
	ngx_rbtree_node_t h;
};
_static_assert(sizeof(struct tree_node_s) == 40);
_static_assert(sizeof(struct tree_node_intl_s) == 40);
_static_assert(offsetof(struct tree_node_s, key) == TREE_OBJECT_OFFSET);
_static_assert(offsetof(struct tree_node_intl_s, h.key) == TREE_OBJECT_OFFSET);
_static_assert(tree_get_object(NULL) == (void *)TREE_OBJECT_OFFSET);


/**
 * @fn tree_clean
 */
void tree_clean(
	tree_t *_tree)
{
	struct tree_s *tree = (struct tree_s *)_tree;
	if(tree == NULL) { return; }

	lmm_t *lmm = tree->lmm;
	struct tree_vec_s *v = tree->vroot;
	while(v != NULL) {
		struct tree_vec_s *vnext = v->next;
		lmm_free(lmm, v); v = vnext;
	}
	lmm_free(lmm, tree);
	return;
}

/**
 * @fn tree_init
 */
tree_t *tree_init(
	uint64_t object_size,
	tree_params_t const *params)
{
	struct tree_params_s const default_params = { 0 };
	params = (params == NULL) ? &default_params : params;

	/* malloc mem */
	lmm_t *lmm = (lmm_t *)params->lmm;
	struct tree_s *tree = (struct tree_s *)lmm_malloc(lmm, sizeof(struct tree_s));
	if(tree == NULL) {
		return(NULL);
	}
	memset(tree, 0, sizeof(struct tree_s));

	/* set params */
	tree->lmm = lmm;
	tree->object_size = _roundup(object_size + sizeof(ngx_rbtree_node_t), 16);
	tree->params = *params;

	/* init vector */
	tree->vhead = tree->vroot = (struct tree_vec_s *)lmm_malloc(lmm,
		sizeof(struct tree_vec_s) + TREE_INIT_ELEM_CNT * tree->object_size);
	tree->vhead->next = NULL;			/* dual linked list */
	tree->vhead->prev = NULL;
	tree->vhead->cnt = TREE_INIT_ELEM_CNT;

	tree->v = tree->vhead->v;
	tree->vrem = tree->vhead->cnt;

	/* init free list */
	struct tree_node_intl_s *tail = (struct tree_node_intl_s *)tree->v;
	tail->h.key = (int64_t)NULL;
	debug("vector inited, v(%p), head(%p), root(%p)", tree->v, tree->vhead, tree->vroot);

	/* init tree */
	ngx_rbtree_init(&tree->t, &tree->sentinel, ngx_rbtree_insert_value);
	tree->sentinel.left = (void *)0x01;
	tree->sentinel.right = (void *)0x02;
	tree->sentinel.data = 0xff;
	return((tree_t *)tree);
}

/**
 * @fn tree_flush
 */
void tree_flush(
	tree_t *_tree)
{
	struct tree_s *tree = (struct tree_s *)_tree;
	if(tree == NULL) { return; }

	/* init v */
	tree->vhead = tree->vroot;

	tree->v = tree->vhead->v;
	tree->vrem = tree->vhead->cnt;

	/* init free list */
	struct tree_node_intl_s *tail = (struct tree_node_intl_s *)tree->v;
	tail->h.key = (int64_t)NULL;
	debug("vector inited, v(%p), head(%p), root(%p)", tree->v, tree->vhead, tree->vroot);

	/* init tree */
	ngx_rbtree_init(&tree->t, &tree->sentinel, ngx_rbtree_insert_value);
	tree->sentinel.data = 0xff;
	return;
}

/**
 * @fn tree_create_node
 *
 * @brief create a new node (not inserted in the tree)
 */
tree_node_t *tree_create_node(
	tree_t *_tree)
{
	struct tree_s *tree = (struct tree_s *)_tree;
	struct tree_node_intl_s *tail = (struct tree_node_intl_s *)tree->v;
	struct tree_node_intl_s *node = NULL;

	/* check the recycle list */
	if((struct tree_node_intl_s *)tail->h.key != NULL) {
		/* recycle removed space */
		node = (struct tree_node_intl_s *)tail->h.key;
		debug("node recycled, node(%p)", node);

		/* update root of the freed list */
		tail->h.key = node->h.key;
	} else {
		/* add new node */
		node = tail;
		tree->v += tree->object_size;

		if(--tree->vrem <= 0) {
			if(tree->vhead->next == NULL) {
				/* add new vector */
				int64_t next_cnt = tree->vhead->cnt * 2;
				struct tree_vec_s *v = tree->vhead->next = (struct tree_vec_s *)lmm_malloc(
					tree->lmm, sizeof(struct tree_vec_s) + next_cnt * tree->object_size);
				v->next = NULL;
				v->prev = tree->vhead;
				v->cnt = next_cnt;
			}

			/* follow the forward link */
			tree->vhead = tree->vhead->next;

			/* init v and vrem */
			tree->v = tree->vhead->v;
			tree->vrem = tree->vhead->cnt;
			debug("added new vector, v(%p), vhead(%p), vroot(%p)", tree->v, tree->vhead, tree->vroot);
		}

		/* copy root of free list */
		struct tree_node_intl_s *tail = (struct tree_node_intl_s *)tree->v;
		tail->h.key = node->h.key;
		debug("new node created, node(%p)", node);
	}

	/* mark node */
	node->h.data = 0xff;
	return((tree_node_t *)node);
}

/**
 * @fn tree_insert
 *
 * @brief insert a node
 */
void tree_insert(
	tree_t *_tree,
	tree_node_t *_node)
{
	struct tree_s *tree = (struct tree_s *)_tree;
	struct tree_node_intl_s *node = (struct tree_node_intl_s *)_node;
	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);
	ngx_rbtree_insert(&tree->t, (ngx_rbtree_node_t *)node);
	return;
}

/**
 * @fn tree_remove
 *
 * @brief remove a node, automatically freed if malloc'd with tree_reserve_node
 */
void tree_remove(
	tree_t *_tree,
	tree_node_t *_node)
{
	struct tree_s *tree = (struct tree_s *)_tree;
	struct tree_node_intl_s *node = (struct tree_node_intl_s *)_node;
	ngx_rbtree_delete(&tree->t, (ngx_rbtree_node_t *)node);

	if(node->h.data == 0xff) {
		/* append node to the head of freed list */
		struct tree_node_intl_s *tail = (struct tree_node_intl_s *)tree->v;
		node->h.key = tail->h.key;
		tail->h.key = (int64_t)node;
	}
	return;
}

/**
 * @fn tree_search_key
 *
 * @brief search a node by key, returning the leftmost node
 */
tree_node_t *tree_search_key(
	tree_t *_tree,
	int64_t key)
{
	struct tree_s *tree = (struct tree_s *)_tree;
	return((tree_node_t *)ngx_rbtree_find_key(&tree->t, key));
}

/**
 * @fn tree_search_key_left
 *
 * @brief search a node by key. returns the nearest node in the left half of the tree if key was not found.
 */
tree_node_t *tree_search_key_left(
	tree_t *_tree,
	int64_t key)
{
	struct tree_s *tree = (struct tree_s *)_tree;
	return((tree_node_t *)ngx_rbtree_find_key_left(&tree->t, key));
}

/**
 * @fn tree_search_key_right
 *
 * @brief search a node by key. returns the nearest node in the right half of the tree if key was not found.
 */
tree_node_t *tree_search_key_right(
	tree_t *_tree,
	int64_t key)
{
	struct tree_s *tree = (struct tree_s *)_tree;
	return((tree_node_t *)ngx_rbtree_find_key_right(&tree->t, key));
}

/**
 * @fn tree_left
 *
 * @brief returns the left next node
 */
tree_node_t *tree_left(
	tree_t *_tree,
	tree_node_t const *node)
{
	struct tree_s *tree = (struct tree_s *)_tree;
	return((tree_node_t *)ngx_rbtree_find_left(&tree->t, (ngx_rbtree_node_t *)node));
}

/**
 * @fn tree_right
 *
 * @brief returns the right next node
 */
tree_node_t *tree_right(
	tree_t *_tree,
	tree_node_t const *node)
{
	struct tree_s *tree = (struct tree_s *)_tree;
	return((tree_node_t *)ngx_rbtree_find_right(&tree->t, (ngx_rbtree_node_t *)node));
}

/**
 * @fn tree_walk
 *
 * @brief walk over the tree
 */
void tree_walk_intl(
	ngx_rbtree_node_t **node,
	ngx_rbtree_node_t *sentinel,
	void *_ctx)
{
	tree_t *tree = (tree_t *)_ctx;
	tree->fn((tree_node_t *)(*node), tree->ctx);
	return;
}
void tree_walk(
	tree_t *_tree,
	tree_walk_t _fn,
	void *_ctx)
{
	struct tree_s *tree = (struct tree_s *)_tree;

	tree->fn = _fn;
	tree->ctx = _ctx;
	ngx_rbtree_walk(&tree->t, (ngx_rbtree_walk_pt)tree_walk_intl, (void *)tree);
	return;
}


/* unittests */
unittest_config(
	.name = "tree"
);

/* create tree object */
unittest()
{
	tree_t *tree = tree_init(8, NULL);

	assert(tree != NULL);

	tree_clean(tree);
}

/* create node */
unittest()
{
	tree_t *tree = tree_init(8, NULL);

	tree_node_t *node = tree_create_node(tree);
	assert(node != NULL);

	int64_t *obj = (int64_t *)tree_get_object(node);
	obj[0] = 0xcafebabe;
	obj[1] = 0x12345678;

	/* insert and search */
	tree_insert(tree, node);
	tree_node_t *found = tree_search_key(tree, 0xcafebabe);
	assert(found == node, "found(%p), node(%p)", found, node);
	assert(obj[0] == 0xcafebabe, "obj[0](%lld)", obj[0]);
	assert(obj[1] == 0x12345678, "obj[1](%lld)", obj[1]);

	/* remove */
	tree_remove(tree, node);
	found = tree_search_key(tree, 0xcafebabe);
	assert(found == NULL, "found(%p), node(%p)", found, node);

	tree_clean(tree);
}

/* create multiple nodes */
unittest()
{
	tree_t *tree = tree_init(8, NULL);

	#define _shuf(x)		( (0xff & ((i<<4) | (i>>4))) )

	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);

	/* insert */
	for(int64_t i = 0; i < 256; i++) {
		tree_node_t *n = tree_create_node(tree);
		assert(n != NULL);

		int64_t *o = (int64_t *)tree_get_object(n);
		o[0] = _shuf(i)<<1;
		o[1] = i;

		tree_insert(tree, n);
	}

	/* search (1) */
	for(int64_t i = 0; i < 256; i++) {
		tree_node_t *n = tree_search_key(tree, i<<1);
		assert(n != NULL);

		int64_t *obj = (int64_t *)tree_get_object(n);
		assert(obj[0] == i<<1, "obj[0](%lld), key(%lld)", obj[0], i<<1);
		assert(obj[1] == _shuf(i), "obj[1](%lld), val(%lld)", obj[1], _shuf(i));
	}

	/* remove */
	for(int64_t i = 0; i < 128; i++) {
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);
		tree_remove(tree, n);
	}

	/* search (2) */
	for(int64_t i = 0; i < 128; i++) {
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n == NULL);
	}
	for(int64_t i = 128; i < 256; i++) {
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		int64_t *obj = (int64_t *)tree_get_object(n);
		assert(obj[0] == _shuf(i)<<1, "obj[0](%lld), key(%lld)", obj[0], _shuf(i)<<1);
		assert(obj[1] == i, "obj[1](%lld), val(%lld)", obj[1], i);
	}

	/* add again */
	for(int64_t i = 0; i < 128; i++) {
		tree_node_t *n = tree_create_node(tree);
		assert(n != NULL);

		int64_t *o = (int64_t *)tree_get_object(n);
		o[0] = _shuf(i)<<1;
		o[1] = 1024 - i;

		tree_insert(tree, n);
	}

	/* search (3) */
	for(int64_t i = 0; i < 128; i++) {
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		int64_t *obj = (int64_t *)tree_get_object(n);
		assert(obj[0] == _shuf(i)<<1, "obj[0](%lld), key(%lld)", obj[0], _shuf(i)<<1);
		assert(obj[1] == 1024 - i, "obj[1](%lld), val(%lld)", obj[1], 1024 - i);
	}
	for(int64_t i = 128; i < 256; i++) {
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		int64_t *obj = (int64_t *)tree_get_object(n);
		assert(obj[0] == _shuf(i)<<1, "obj[0](%lld), key(%lld)", obj[0], _shuf(i)<<1);
		assert(obj[1] == i, "obj[1](%lld), val(%lld)", obj[1], i);
	}

	tree_clean(tree);
}

/* create multiple nodes with malloc */
void unittest_free_node(tree_node_t *node, void *ctx)
{
	free(node);
	return;
}
unittest()
{
	tree_t *tree = tree_init(8, NULL);

	#define _shuf(x)		( (0xff & ((i<<4) | (i>>4))) )

	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);

	/* insert */
	for(int64_t i = 0; i < 256; i++) {
		tree_node_t *n = (tree_node_t *)malloc(sizeof(tree_node_t) + 8);
		assert(n != NULL);

		memset(n, 0, sizeof(tree_node_t));
		int64_t *o = (int64_t *)tree_get_object(n);
		o[0] = _shuf(i)<<1;
		o[1] = i;

		tree_insert(tree, n);
	}

	/* search (1) */
	for(int64_t i = 0; i < 256; i++) {
		tree_node_t *n = tree_search_key(tree, i<<1);
		assert(n != NULL);

		int64_t *obj = (int64_t *)tree_get_object(n);
		assert(obj[0] == i<<1, "obj[0](%lld), key(%lld)", obj[0], i<<1);
		assert(obj[1] == _shuf(i), "obj[1](%lld), val(%lld)", obj[1], _shuf(i));
	}

	/* remove */
	for(int64_t i = 0; i < 128; i++) {
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);
		tree_remove(tree, n);
		free(n);
	}

	/* search (2) */
	for(int64_t i = 0; i < 128; i++) {
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n == NULL);
	}
	for(int64_t i = 128; i < 256; i++) {
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		int64_t *obj = (int64_t *)tree_get_object(n);
		assert(obj[0] == _shuf(i)<<1, "obj[0](%lld), key(%lld)", obj[0], _shuf(i)<<1);
		assert(obj[1] == i, "obj[1](%lld), val(%lld)", obj[1], i);
	}

	/* add again */
	for(int64_t i = 0; i < 128; i++) {
		tree_node_t *n = (tree_node_t *)malloc(sizeof(tree_node_t) + 8);
		assert(n != NULL);

		memset(n, 0, sizeof(tree_node_t));
		int64_t *o = (int64_t *)tree_get_object(n);
		o[0] = _shuf(i)<<1;
		o[1] = 1024 - i;

		tree_insert(tree, n);
	}

	/* search (3) */
	for(int64_t i = 0; i < 128; i++) {
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		int64_t *obj = (int64_t *)tree_get_object(n);
		assert(obj[0] == _shuf(i)<<1, "obj[0](%lld), key(%lld)", obj[0], _shuf(i)<<1);
		assert(obj[1] == 1024 - i, "obj[1](%lld), val(%lld)", obj[1], 1024 - i);
	}
	for(int64_t i = 128; i < 256; i++) {
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		int64_t *obj = (int64_t *)tree_get_object(n);
		assert(obj[0] == _shuf(i)<<1, "obj[0](%lld), key(%lld)", obj[0], _shuf(i)<<1);
		assert(obj[1] == i, "obj[1](%lld), val(%lld)", obj[1], i);
	}

	/* cleanup */
	tree_walk(tree, unittest_free_node, NULL);
	tree_clean(tree);
}

/* flush, duplicated key */
unittest()
{
	tree_t *tree = tree_init(8, NULL);

	#define _shuf(x)		( (0xff & ((i<<4) | (i>>4))) )

	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);

	/* insert */
	for(int64_t i = 0; i < 65536; i++) {
		debug("insert node i(%lld)", i);
		tree_node_t *n = tree_create_node(tree);
		assert(n != NULL);

		int64_t *o = (int64_t *)tree_get_object(n);
		o[0] = _shuf(i)<<1;
		o[1] = 65536 - i;

		tree_insert(tree, n);
	}

	/* flush and insert again */
	tree_flush(tree);
	for(int64_t i = 0; i < 65536; i++) {
		debug("insert node again i(%lld)", i);
		tree_node_t *n = tree_create_node(tree);
		assert(n != NULL);

		int64_t *o = (int64_t *)tree_get_object(n);
		o[0] = _shuf(i)<<1;
		o[1] = i;

		tree_insert(tree, n);
	}

	/* search (1) */
	for(int64_t i = 0; i < 256; i++) {
		tree_node_t *n = tree_search_key(tree, i<<1);
		assert(n != NULL);
	}
	for(int64_t i = 256; i < 65536; i++) {
		tree_node_t *n = tree_search_key(tree, i<<1);
		assert(n == NULL);
	}

	/* remove */
	for(int64_t i = 0; i < 65536; i++) {
		debug("remove node i(%lld)", i);
		tree_node_t *n = tree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);
		tree_remove(tree, n);
	}

	/* search (2) */
	for(int64_t i = 0; i < 256; i++) {
		tree_node_t *n = tree_search_key(tree, i<<1);
		assert(n == NULL);
	}
	tree_clean(tree);
}

/* search_key_left / search_key_right / left / right (NULL) */
unittest()
{
	tree_t *tree = tree_init(8, NULL);

	#define _shuf(x)		( (0xff & ((i<<4) | (i>>4))) )

	/* search left */
	tree_node_t *n = tree_search_key_left(tree, 0);
	assert(n == NULL);

	n = tree_search_key_left(tree, 100);
	assert(n == NULL);

	/* search right */
	n = tree_search_key_right(tree, 0);
	assert(n == NULL);

	n = tree_search_key_right(tree, 200);
	assert(n == NULL);

	/* left / right of sentinel */
	n = tree_left(tree, (tree_node_t *)&tree->sentinel);
	assert(n == NULL);

	n = tree_right(tree, (tree_node_t *)&tree->sentinel);
	assert(n == NULL);

	tree_clean(tree);
}

/* search_key_left / search_key_right / left / right (NULL) */
unittest()
{
	tree_t *tree = tree_init(8, NULL);

	#define _shuf(x)		( (0xff & ((i<<4) | (i>>4))) )

	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);

	/* insert */
	for(int64_t i = 0; i < 256; i++) {
		tree_node_t *n = tree_create_node(tree);
		assert(n != NULL);

		int64_t *o = (int64_t *)tree_get_object(n);
		o[0] = _shuf(i)<<1;
		o[1] = i;

		tree_insert(tree, n);
	}

	/* search left */
	tree_node_t *n1 = tree_search_key_left(tree, 64);
	tree_node_t *n2 = tree_search_key_left(tree, 65);
	assert(n1 == n2, "n1(%p), n2(%p)", n1, n2);
	assert(n1->key == 64, "key(%lld)", n1->key);
	assert(n2->key == 64, "key(%lld)", n2->key);

	/* left and right */
	n1 = tree_left(tree, n1);
	assert(n1->key == 62, "key(%lld)", n1->key);
	n2 = tree_right(tree, n2);
	assert(n2->key == 66, "key(%lld)", n2->key);


	/* search right */
	n1 = tree_search_key_right(tree, 31);
	n2 = tree_search_key_right(tree, 32);
	assert(n1 == n2, "n1(%p), n2(%p)", n1, n2);
	assert(n1->key == 32, "key(%lld)", n1->key);
	assert(n2->key == 32, "key(%lld)", n2->key);

	/* left and right */
	n1 = tree_left(tree, n1);
	assert(n1->key == 30, "key(%lld)", n1->key);
	n2 = tree_right(tree, n2);
	assert(n2->key == 34, "key(%lld)", n2->key);

	tree_clean(tree);
}

/**
 * end of tree.c
 */
