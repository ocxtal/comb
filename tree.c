
/**
 * @file tree.c
 *
 * @brief an wrapper of ngx_rbtree.c in nginx (https://nginx.org/) core library
 */

#define UNITTEST_UNIQUE_ID		59
#include "unittest.h"

#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include "ngx_rbtree.h"
#include "lmm.h"
#include "log.h"
#include "sassert.h"
#include "tree.h"


/* constants */
#define	RBTREE_INIT_ELEM_CNT		( 64 )

/* roundup */
#define _roundup(x, base)			( ((x) + (base) - 1) & ~((base) - 1) )

/**
 * @struct rbtree_vec_s
 */
struct rbtree_vec_s {
	struct rbtree_vec_s *next;
	struct rbtree_vec_s *prev;
	int64_t cnt;
	int64_t pad;
	uint8_t v[];
};

/**
 * @struct rbtree_s
 */
struct rbtree_s {
	lmm_t *lmm;
	lmm_t *lmm_iter;
	uint32_t object_size;
	uint32_t pad;
	struct rbtree_params_s params;

	/* vector pointers */
	uint8_t *v;
	int64_t vrem;
	struct rbtree_vec_s *vhead, *vroot;

	/* working buffer */
	rbtree_walk_t fn;
	void *ctx;

	/* tree */
	ngx_rbtree_t t;
	ngx_rbtree_node_t sentinel;

	/* reserved */
	uint8_t reserved[sizeof(struct ngx_ivtree_node_s)
		- sizeof(ngx_rbtree_node_t)];
};

/**
 * @struct ivtree_iter_s
 */
struct ivtree_iter_s {
	lmm_t *lmm;
	ngx_rbtree_t *t;
	int64_t llim, rlim, tlim;
	ngx_ivtree_node_t *node;
};


/* assertions */
_static_assert(sizeof(struct rbtree_node_s) == 40);
_static_assert(sizeof(ngx_rbtree_node_t) == 40);


/**
 * @fn rbtree_clean
 */
void rbtree_clean(
	rbtree_t *_tree)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	if(tree == NULL) { return; }

	lmm_t *lmm = tree->lmm;
	struct rbtree_vec_s *v = tree->vroot;
	while(v != NULL) {
		struct rbtree_vec_s *vnext = v->next;
		lmm_free(lmm, v); v = vnext;
	}
	lmm_free(lmm, tree->lmm_iter);
	lmm_free(lmm, tree);
	return;
}

/**
 * @fn rbtree_init
 */
rbtree_t *rbtree_init(
	uint64_t object_size,
	rbtree_params_t const *params)
{
	struct rbtree_params_s const default_params = { 0 };
	params = (params == NULL) ? &default_params : params;

	/* malloc mem */
	lmm_t *lmm = (lmm_t *)params->lmm;
	struct rbtree_s *tree = (struct rbtree_s *)lmm_malloc(lmm, sizeof(struct rbtree_s));
	if(tree == NULL) {
		return(NULL);
	}
	memset(tree, 0, sizeof(struct rbtree_s));
	lmm_t *lmm_iter = lmm_init(
		lmm_malloc(lmm, 3 * sizeof(struct ivtree_iter_s)),
		3 * sizeof(struct ivtree_iter_s));

	/* set params */
	tree->lmm = lmm;
	tree->lmm_iter = lmm_iter;
	tree->object_size = _roundup(object_size, 16);
	tree->params = *params;

	/* init vector */
	tree->vhead = tree->vroot = (struct rbtree_vec_s *)lmm_malloc(lmm,
		sizeof(struct rbtree_vec_s) + RBTREE_INIT_ELEM_CNT * tree->object_size);
	tree->vhead->next = NULL;			/* dual linked list */
	tree->vhead->prev = NULL;
	tree->vhead->cnt = RBTREE_INIT_ELEM_CNT;

	tree->v = tree->vhead->v;
	tree->vrem = tree->vhead->cnt;

	/* init free list */
	ngx_rbtree_node_t *tail = (ngx_rbtree_node_t *)tree->v;
	tail->key = (int64_t)NULL;
	debug("vector inited, v(%p), head(%p), root(%p)", tree->v, tree->vhead, tree->vroot);

	/* init tree */
	ngx_rbtree_init(&tree->t, &tree->sentinel, ngx_rbtree_insert_value);
	tree->sentinel.left = (void *)0x01;
	tree->sentinel.right = (void *)0x02;
	tree->sentinel.data = 0xff;
	return((rbtree_t *)tree);
}

/**
 * @fn rbtree_flush
 */
void rbtree_flush(
	rbtree_t *_tree)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	if(tree == NULL) { return; }

	/* init v */
	tree->vhead = tree->vroot;

	tree->v = tree->vhead->v;
	tree->vrem = tree->vhead->cnt;

	/* init free list */
	ngx_rbtree_node_t *tail = (ngx_rbtree_node_t *)tree->v;
	tail->key = (int64_t)NULL;
	debug("vector inited, v(%p), head(%p), root(%p)", tree->v, tree->vhead, tree->vroot);

	/* init tree */
	ngx_rbtree_init(&tree->t, &tree->sentinel, ngx_rbtree_insert_value);
	tree->sentinel.data = 0xff;
	return;
}

/**
 * @fn rbtree_create_node
 *
 * @brief create a new node (not inserted in the tree)
 */
rbtree_node_t *rbtree_create_node(
	rbtree_t *_tree)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	ngx_rbtree_node_t *tail = (ngx_rbtree_node_t *)tree->v;
	ngx_rbtree_node_t *node = NULL;

	/* check the recycle list */
	if((ngx_rbtree_node_t *)tail->key != NULL) {
		/* recycle removed space */
		node = (ngx_rbtree_node_t *)tail->key;
		debug("node recycled, node(%p)", node);

		/* update root of the freed list */
		tail->key = node->key;
	} else {
		/* add new node */
		node = tail;
		tree->v += tree->object_size;

		if(--tree->vrem <= 0) {
			if(tree->vhead->next == NULL) {
				/* add new vector */
				int64_t next_cnt = tree->vhead->cnt * 2;
				struct rbtree_vec_s *v = tree->vhead->next = (struct rbtree_vec_s *)lmm_malloc(
					tree->lmm, sizeof(struct rbtree_vec_s) + next_cnt * tree->object_size);
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
		ngx_rbtree_node_t *tail = (ngx_rbtree_node_t *)tree->v;
		tail->key = node->key;
		debug("new node created, node(%p)", node);
	}

	/* mark node */
	node->data = 0xff;
	return((rbtree_node_t *)node);
}

/**
 * @fn rbtree_insert
 *
 * @brief insert a node
 */
void rbtree_insert(
	rbtree_t *_tree,
	rbtree_node_t *_node)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	ngx_rbtree_node_t *node = (ngx_rbtree_node_t *)_node;
	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);
	ngx_rbtree_insert(&tree->t, node);
	return;
}

/**
 * @fn rbtree_remove
 *
 * @brief remove a node, automatically freed if malloc'd with rbtree_reserve_node
 */
void rbtree_remove(
	rbtree_t *_tree,
	rbtree_node_t *_node)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	ngx_rbtree_node_t *node = (ngx_rbtree_node_t *)_node;
	ngx_rbtree_delete(&tree->t, node);

	if(node->data == 0xff) {
		/* append node to the head of freed list */
		ngx_rbtree_node_t *tail = (ngx_rbtree_node_t *)tree->v;
		node->key = tail->key;
		tail->key = (int64_t)node;
	}
	return;
}

/**
 * @fn rbtree_search_key
 *
 * @brief search a node by key, returning the leftmost node
 */
rbtree_node_t *rbtree_search_key(
	rbtree_t *_tree,
	int64_t key)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	return((rbtree_node_t *)ngx_rbtree_find_key(&tree->t, key));
}

/**
 * @fn rbtree_search_key_left
 *
 * @brief search a node by key. returns the nearest node in the left half of the tree if key was not found.
 */
rbtree_node_t *rbtree_search_key_left(
	rbtree_t *_tree,
	int64_t key)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	return((rbtree_node_t *)ngx_rbtree_find_key_left(&tree->t, key));
}

/**
 * @fn rbtree_search_key_right
 *
 * @brief search a node by key. returns the nearest node in the right half of the tree if key was not found.
 */
rbtree_node_t *rbtree_search_key_right(
	rbtree_t *_tree,
	int64_t key)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	return((rbtree_node_t *)ngx_rbtree_find_key_right(&tree->t, key));
}

/**
 * @fn rbtree_left
 *
 * @brief returns the left next node
 */
rbtree_node_t *rbtree_left(
	rbtree_t *_tree,
	rbtree_node_t const *node)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	return((rbtree_node_t *)ngx_rbtree_find_left(&tree->t, (ngx_rbtree_node_t *)node));
}

/**
 * @fn rbtree_right
 *
 * @brief returns the right next node
 */
rbtree_node_t *rbtree_right(
	rbtree_t *_tree,
	rbtree_node_t const *node)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	return((rbtree_node_t *)ngx_rbtree_find_right(&tree->t, (ngx_rbtree_node_t *)node));
}

/**
 * @fn rbtree_walk
 *
 * @brief walk over the tree
 */
void rbtree_walk_intl(
	ngx_rbtree_node_t **node,
	ngx_rbtree_node_t *sentinel,
	void *_ctx)
{
	rbtree_t *tree = (rbtree_t *)_ctx;
	tree->fn((rbtree_node_t *)(*node), tree->ctx);
	return;
}
void rbtree_walk(
	rbtree_t *_tree,
	rbtree_walk_t _fn,
	void *_ctx)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;

	tree->fn = _fn;
	tree->ctx = _ctx;
	ngx_rbtree_walk(&tree->t, (ngx_rbtree_walk_pt)rbtree_walk_intl, (void *)tree);
	return;
}


/* interval tree implementation */
/**
 * @fn ivtree_clean
 */
void ivtree_clean(
	ivtree_t *tree)
{
	rbtree_clean((struct rbtree_s *)tree);
	return;
}

/**
 * @fn ivtree_init
 */
ivtree_t *ivtree_init(
	uint64_t object_size,
	ivtree_params_t const *params)
{
	struct rbtree_s *tree = rbtree_init(object_size, (rbtree_params_t const *)params);
	struct ngx_ivtree_node_s *sentinel = (struct ngx_ivtree_node_s *)(&tree->sentinel);

	sentinel->lkey = INT64_MIN;
	sentinel->rkey = INT64_MIN;
	sentinel->rkey_max = INT64_MIN;
	return((ivtree_t *)tree);
}

/**
 * @fn ivtree_flush
 */
void ivtree_flush(
	ivtree_t *tree)
{
	rbtree_flush((rbtree_t *)tree);
	return;
}

/**
 * @fn ivtree_create_node
 *
 * @brief create a new node (not inserted in the tree)
 */
ivtree_node_t *ivtree_create_node(
	ivtree_t *tree)
{
	return((ivtree_node_t *)rbtree_create_node((rbtree_t *)tree));
}

/**
 * @fn ivtree_insert
 *
 * @brief insert a node
 */
void ivtree_insert(
	ivtree_t *_tree,
	ivtree_node_t *_node)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	struct ngx_ivtree_node_s *node = (struct ngx_ivtree_node_s *)_node;
	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);
	ngx_ivtree_insert((ngx_ivtree_t *)&tree->t, node);
	return;
}

/**
 * @fn ivtree_remove
 *
 * @brief remove a node, automatically freed if malloc'd with ivtree_reserve_node
 */
void ivtree_remove(
	ivtree_t *_tree,
	ivtree_node_t *_node)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	ngx_rbtree_node_t *node = (ngx_rbtree_node_t *)_node;
	ngx_ivtree_delete(
		(ngx_ivtree_t *)&tree->t,
		(ngx_ivtree_node_t *)node);

	if(node->data == 0xff) {
		/* append node to the head of freed list */
		ngx_rbtree_node_t *tail = (ngx_rbtree_node_t *)tree->v;
		node->key = tail->key;
		tail->key = (int64_t)node;
	}
	return;
}

/**
 * @fn ivtree_next_node
 */
static inline
ngx_ivtree_node_t *ivtree_next_node(
	ngx_rbtree_t *t,
	ngx_ivtree_node_t *node,
	int64_t llim,
	int64_t rlim,
	int64_t tlim)
{
	debug("ivtree_next_node, llim(%lld), rlim(%lld), tlim(%lld)", llim, rlim, tlim);

	while(node != NULL) {
        debug("check node(%p, %lld, %lld, %lld)",
            node, node->lkey, node->rkey, node->rkey_max);

		if(node->lkey >= tlim) { node = NULL; break; }
		if((uint64_t)(node->rkey - llim) < (rlim - llim)) { break; }

		/* not found, get next */
		node = (ngx_ivtree_node_t *)ngx_rbtree_find_right(
			t, (ngx_rbtree_node_t *)node);
	}
	return(node);
}

/**
 * @fn ivtree_contained
 * @brief return a set of sections contained in [lkey, rkey)
 */
ivtree_iter_t *ivtree_contained(
	ivtree_t *_tree,
	int64_t lkey,
	int64_t rkey)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	struct ivtree_iter_s *iter = (struct ivtree_iter_s *)lmm_malloc(
		tree->lmm_iter, sizeof(struct ivtree_iter_s));
	*iter = (struct ivtree_iter_s){
		.t = &tree->t,
		.lmm = tree->lmm_iter,
		.llim = INT64_MIN,
		.rlim = rkey,
		.tlim = rkey,
		.node = (ngx_ivtree_node_t *)ngx_rbtree_find_key_right(&tree->t, lkey)
	};
	return(iter);
}

/**
 * @fn ivtree_containing
 * @brief return a set of sections containing [lkey, rkey)
 */
ivtree_iter_t *ivtree_containing(
	ivtree_t *_tree,
	int64_t lkey,
	int64_t rkey)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	struct ivtree_iter_s *iter = (struct ivtree_iter_s *)lmm_malloc(
		tree->lmm_iter, sizeof(struct ivtree_iter_s));

	ngx_ivtree_node_t *node = (ngx_ivtree_node_t *)tree->t.root;
	ngx_ivtree_node_t *sentinel = (ngx_ivtree_node_t *)tree->t.sentinel;
	while(node != sentinel && node->left->rkey_max >= rkey) {
		node = node->left;
	}
	*iter = (struct ivtree_iter_s){
		.t = &tree->t,
		.lmm = tree->lmm_iter,
		.llim = rkey,
		.rlim = INT64_MAX,
		.tlim = lkey + 1,
		.node = node
	};
	return(iter);
}

/**
 * @fn ivtree_intersect
 * @brief return a set of sections intersect with [lkey, rkey)
 */
ivtree_iter_t *ivtree_intersect(
	ivtree_t *_tree,
	int64_t lkey,
	int64_t rkey)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;
	struct ivtree_iter_s *iter = (struct ivtree_iter_s *)lmm_malloc(
		tree->lmm_iter, sizeof(struct ivtree_iter_s));

	ngx_ivtree_node_t *node = (ngx_ivtree_node_t *)tree->t.root;
	ngx_ivtree_node_t *sentinel = (ngx_ivtree_node_t *)tree->t.sentinel;
	while(node != sentinel && node->left->rkey_max > lkey) {
		node = node->left;
	}
	*iter = (struct ivtree_iter_s){
		.t = &tree->t,
		.lmm = tree->lmm_iter,
		.llim = lkey + 1,
		.rlim = INT64_MAX,
		.tlim = rkey,
		.node = node
	};
	return(iter);
}

/**
 * @fn ivtree_next
 */
ivtree_node_t *ivtree_next(
	ivtree_iter_t *_iter)
{
	struct ivtree_iter_s *iter = (struct ivtree_iter_s *)_iter;
	ngx_ivtree_node_t *node = ivtree_next_node(iter->t,
		iter->node, iter->llim, iter->rlim, iter->tlim);
	if(node == NULL) { return(NULL); }

	iter->node = (ngx_ivtree_node_t *)ngx_rbtree_find_right(
		iter->t, (ngx_rbtree_node_t *)node);
	return((ivtree_node_t *)node);
}

/**
 * @fn ivtree_iter_clean
 */
void ivtree_iter_clean(
	ivtree_iter_t *_iter)
{
	struct ivtree_iter_s *iter = (struct ivtree_iter_s *)_iter;
	if(iter == NULL) { return; }
	lmm_t *lmm_iter = iter->lmm;
	lmm_free(lmm_iter, iter);
	return;
}

/**
 * @fn ivtree_walk
 *
 * @brief walk over the tree
 */
void ivtree_walk(
	ivtree_t *_tree,
	ivtree_walk_t _fn,
	void *_ctx)
{
	struct rbtree_s *tree = (struct rbtree_s *)_tree;

	tree->fn = (rbtree_walk_t)_fn;
	tree->ctx = _ctx;
	ngx_rbtree_walk(&tree->t, (ngx_rbtree_walk_pt)rbtree_walk_intl, (void *)tree);
	return;
}


/* unittests */
unittest_config(
	.name = "tree"
);

/**
 * @struct ut_rbnode_s
 */
struct ut_rbnode_s {
	rbtree_node_t h;
	int64_t val;
};

/* create tree object */
unittest()
{
	rbtree_t *tree = rbtree_init(sizeof(struct ut_rbnode_s), NULL);

	assert(tree != NULL);

	rbtree_clean(tree);
}

/* create node */
unittest()
{
	rbtree_t *tree = rbtree_init(sizeof(struct ut_rbnode_s), NULL);
	struct ut_rbnode_s *n = (struct ut_rbnode_s *)
		rbtree_create_node(tree);
	assert(n != NULL);

	n->h.key = 0xcafebabe;
	n->val = 0x12345678;

	/* insert and search */
	rbtree_insert(tree, (rbtree_node_t *)n);
	struct ut_rbnode_s *found = (struct ut_rbnode_s *)
		 rbtree_search_key(tree, 0xcafebabe);
	assert(found == n, "found(%p), n(%p)", found, n);
	assert(n->h.key == 0xcafebabe, "n->h.key(%lld)", n->h.key);
	assert(n->val == 0x12345678, "n->val(%lld)", n->val);

	/* remove */
	rbtree_remove(tree, (rbtree_node_t *)n);
	found = (struct ut_rbnode_s *)rbtree_search_key(tree, 0xcafebabe);
	assert(found == NULL, "found(%p), n(%p)", found, n);

	rbtree_clean(tree);
}

/* create multiple nodes */
unittest()
{
	rbtree_t *tree = rbtree_init(sizeof(struct ut_rbnode_s), NULL);

	#define _shuf(x)		( (0xff & ((i<<4) | (i>>4))) )

	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);

	/* insert */
	for(int64_t i = 0; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_create_node(tree);
		assert(n != NULL);

		n->h.key = _shuf(i)<<1;
		n->val = i;
		rbtree_insert(tree, (rbtree_node_t *)n);
	}

	/* search (1) */
	for(int64_t i = 0; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, i<<1);
		assert(n != NULL);

		assert(n->h.key == i<<1, "n->h.key(%lld), key(%lld)", n->h.key, i<<1);
		assert(n->val == _shuf(i), "n->val(%lld), val(%lld)", n->val, _shuf(i));
	}

	/* remove */
	for(int64_t i = 0; i < 128; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);
		rbtree_remove(tree, (rbtree_node_t *)n);
	}

	/* search (2) */
	for(int64_t i = 0; i < 128; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n == NULL);
	}
	for(int64_t i = 128; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		assert(n->h.key == _shuf(i)<<1, "n->h.key(%lld), key(%lld)", n->h.key, _shuf(i)<<1);
		assert(n->val == i, "n->val(%lld), val(%lld)", n->val, i);
	}

	/* add again */
	for(int64_t i = 0; i < 128; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_create_node(tree);
		assert(n != NULL);

		n->h.key = _shuf(i)<<1;
		n->val = 1024 - i;
		rbtree_insert(tree, (rbtree_node_t *)n);
	}

	/* search (3) */
	for(int64_t i = 0; i < 128; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		assert(n->h.key == _shuf(i)<<1, "n->h.key(%lld), key(%lld)", n->h.key, _shuf(i)<<1);
		assert(n->val == 1024 - i, "n->val(%lld), val(%lld)", n->val, 1024 - i);
	}
	for(int64_t i = 128; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		assert(n->h.key == _shuf(i)<<1, "n->h.key(%lld), key(%lld)", n->h.key, _shuf(i)<<1);
		assert(n->val == i, "n->val(%lld), val(%lld)", n->val, i);
	}

	rbtree_clean(tree);
}

unittest()
{
	int64_t const cnt = 32 * 1024 * 1024;
	srand(time(NULL));

	rbtree_t *tree = rbtree_init(sizeof(struct ut_rbnode_s), NULL);

	/* insert */
	for(int64_t i = 0; i < cnt; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_create_node(tree);
		assert(n != NULL);

		n->h.key = rand();
		n->val = i;
		rbtree_insert(tree, (rbtree_node_t *)n);
	}

	/* remove */
	for(int64_t i = 0; i < cnt; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, i);
		if(n != NULL) {
			rbtree_remove(tree, (rbtree_node_t *)n);
		}
	}
	rbtree_clean(tree);
}

/* create multiple nodes with malloc */
void unittest_free_node(rbtree_node_t *node, void *ctx)
{
	free(node);
	return;
}
unittest()
{
	rbtree_t *tree = rbtree_init(sizeof(struct ut_rbnode_s), NULL);

	#define _shuf(x)		( (0xff & ((i<<4) | (i>>4))) )

	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);

	/* insert */
	for(int64_t i = 0; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			(rbtree_node_t *)malloc(sizeof(rbtree_node_t) + 8);
		assert(n != NULL);

		memset(n, 0, sizeof(rbtree_node_t));
		n->h.key = _shuf(i)<<1;
		n->val = i;
		rbtree_insert(tree, (rbtree_node_t *)n);
	}

	/* search (1) */
	for(int64_t i = 0; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, i<<1);
		assert(n != NULL);

		assert(n->h.key == i<<1, "n->h.key(%lld), key(%lld)", n->h.key, i<<1);
		assert(n->val == _shuf(i), "n->val(%lld), val(%lld)", n->val, _shuf(i));
	}

	/* remove */
	for(int64_t i = 0; i < 128; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);
		rbtree_remove(tree, (rbtree_node_t *)n);
		free(n);
	}

	/* search (2) */
	for(int64_t i = 0; i < 128; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n == NULL);
	}
	for(int64_t i = 128; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		assert(n->h.key == _shuf(i)<<1, "n->h.key(%lld), key(%lld)", n->h.key, _shuf(i)<<1);
		assert(n->val == i, "n->val(%lld), val(%lld)", n->val, i);
	}

	/* add again */
	for(int64_t i = 0; i < 128; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			(rbtree_node_t *)malloc(sizeof(rbtree_node_t) + 8);
		assert(n != NULL);

		memset(n, 0, sizeof(rbtree_node_t));
		n->h.key = _shuf(i)<<1;
		n->val = 1024 - i;
		rbtree_insert(tree, (rbtree_node_t *)n);
	}

	/* search (3) */
	for(int64_t i = 0; i < 128; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		assert(n->h.key == _shuf(i)<<1, "n->h.key(%lld), key(%lld)", n->h.key, _shuf(i)<<1);
		assert(n->val == 1024 - i, "n->val(%lld), val(%lld)", n->val, 1024 - i);
	}
	for(int64_t i = 128; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);

		assert(n->h.key == _shuf(i)<<1, "n->h.key(%lld), key(%lld)", n->h.key, _shuf(i)<<1);
		assert(n->val == i, "n->val(%lld), val(%lld)", n->val, i);
	}

	/* cleanup */
	rbtree_walk(tree, unittest_free_node, NULL);
	rbtree_clean(tree);
}

/* flush, duplicated key */
unittest()
{
	rbtree_t *tree = rbtree_init(sizeof(struct ut_rbnode_s), NULL);

	#define _shuf(x)		( (0xff & ((i<<4) | (i>>4))) )

	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);

	/* insert */
	for(int64_t i = 0; i < 65536; i++) {
		debug("insert node i(%lld)", i);
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_create_node(tree);
		assert(n != NULL);

		n->h.key = _shuf(i)<<1;
		n->val = 65536 - i;
		rbtree_insert(tree, (rbtree_node_t *)n);
	}

	/* flush and insert again */
	rbtree_flush(tree);
	for(int64_t i = 0; i < 65536; i++) {
		debug("insert node again i(%lld)", i);
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_create_node(tree);
		assert(n != NULL);

		n->h.key = _shuf(i)<<1;
		n->val = i;
		rbtree_insert(tree, (rbtree_node_t *)n);
	}

	/* search (1) */
	for(int64_t i = 0; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, i<<1);
		assert(n != NULL);
	}
	for(int64_t i = 256; i < 65536; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, i<<1);
		assert(n == NULL);
	}

	/* remove */
	for(int64_t i = 0; i < 65536; i++) {
		debug("remove node i(%lld)", i);
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, _shuf(i)<<1);
		assert(n != NULL);
		rbtree_remove(tree, (rbtree_node_t *)n);
	}

	/* search (2) */
	for(int64_t i = 0; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_search_key(tree, i<<1);
		assert(n == NULL);
	}
	rbtree_clean(tree);
}

/* search_key_left / search_key_right / left / right (NULL) */
unittest()
{
	rbtree_t *tree = rbtree_init(sizeof(struct ut_rbnode_s), NULL);

	#define _shuf(x)		( (0xff & ((i<<4) | (i>>4))) )

	/* search left */
	struct ut_rbnode_s *n = (struct ut_rbnode_s *)
		rbtree_search_key_left(tree, 0);
	assert(n == NULL);

	n = (struct ut_rbnode_s *)rbtree_search_key_left(tree, 100);
	assert(n == NULL);

	/* search right */
	n = (struct ut_rbnode_s *)rbtree_search_key_right(tree, 0);
	assert(n == NULL);

	n = (struct ut_rbnode_s *)rbtree_search_key_right(tree, 200);
	assert(n == NULL);

	/* left / right of sentinel */
	n = (struct ut_rbnode_s *)rbtree_left(tree, (rbtree_node_t *)&tree->sentinel);
	assert(n == NULL);

	n = (struct ut_rbnode_s *)rbtree_right(tree, (rbtree_node_t *)&tree->sentinel);
	assert(n == NULL);

	rbtree_clean(tree);
}

/* search_key_left / search_key_right / left / right (NULL) */
unittest()
{
	rbtree_t *tree = rbtree_init(sizeof(struct ut_rbnode_s), NULL);

	#define _shuf(x)		( (0xff & ((i<<4) | (i>>4))) )

	debug("tree->root(%p), tree->sentinel(%p)", tree->t.root, tree->t.sentinel);

	/* insert */
	for(int64_t i = 0; i < 256; i++) {
		struct ut_rbnode_s *n = (struct ut_rbnode_s *)
			rbtree_create_node(tree);
		assert(n != NULL);

		n->h.key = _shuf(i)<<1;
		n->val = i;
		rbtree_insert(tree, (rbtree_node_t *)n);
	}

	/* search left */
	rbtree_node_t *n1 = rbtree_search_key_left(tree, 64);
	rbtree_node_t *n2 = rbtree_search_key_left(tree, 65);
	assert(n1 == n2, "n1(%p), n2(%p)", n1, n2);
	assert(n1->key == 64, "key(%lld)", n1->key);
	assert(n2->key == 64, "key(%lld)", n2->key);

	/* left and right */
	n1 = rbtree_left(tree, n1);
	assert(n1->key == 62, "key(%lld)", n1->key);
	n2 = rbtree_right(tree, n2);
	assert(n2->key == 66, "key(%lld)", n2->key);


	/* search right */
	n1 = rbtree_search_key_right(tree, 31);
	n2 = rbtree_search_key_right(tree, 32);
	assert(n1 == n2, "n1(%p), n2(%p)", n1, n2);
	assert(n1->key == 32, "key(%lld)", n1->key);
	assert(n2->key == 32, "key(%lld)", n2->key);

	/* left and right */
	n1 = rbtree_left(tree, n1);
	assert(n1->key == 30, "key(%lld)", n1->key);
	n2 = rbtree_right(tree, n2);
	assert(n2->key == 34, "key(%lld)", n2->key);

	rbtree_clean(tree);
}

/* interval tree test */
/**
 * @struct ut_ivnode_s
 */
struct ut_ivnode_s {
	ivtree_node_t h;
	int64_t val;
};

/* create context */
unittest()
{
	ivtree_t *tree = ivtree_init(sizeof(struct ut_ivnode_s), NULL);
	assert(tree != NULL);
	ivtree_clean(tree);
}

/* insert and remove */
unittest()
{
	ivtree_t *tree = ivtree_init(sizeof(struct ut_ivnode_s), NULL);

	int64_t const cnt = 10;
	int64_t const iv[10][2] = {
		{ 100, 200 },
		{ 50, 100 },
		{ 120, 180 },
		{ 60, 240 },
		{ 70, 200 },
		{ 100, 140 },
		{ 90, 150 },
		{ 110, 160 },
		{ 120, 140 },
		{ 40, 190 }
	};

	/* insert */
	for(int64_t i = 0; i < cnt; i++) {
		struct ut_ivnode_s *n = (struct ut_ivnode_s *)
			ivtree_create_node(tree);
		assert(n != NULL);

		n->h.lkey = iv[i][0];
		n->h.rkey = iv[i][1];
		n->val = 0x12345678;

		ivtree_insert(tree, (ivtree_node_t *)n);
	}

	/* cleanup */
	ivtree_clean(tree);
}

#define _check_node(n, lk, rk) \
	(n) != NULL && (n)->lkey == (lk) && (n)->rkey == (rk)
#define _print_node(n) \
	"node(%p), lkey(%lld), rkey(%lld)", \
	(n), (n)->lkey, (n)->rkey

void ivtree_print_node(
	ivtree_node_t *node,
	void *ctx)
{
	// ngx_ivtree_node_t *n = (ngx_ivtree_node_t *)node;
	// debug("node(%p), left(%p), right(%p), parent(%p), lkey(%lld), rkey(%lld), rkey_max(%lld)",
		// n, n->left, n->right, n->parent, n->lkey, n->rkey, n->rkey_max);
	return;
}

/* contained iterator */
unittest()
{
	ivtree_t *tree = ivtree_init(sizeof(struct ut_ivnode_s), NULL);

	int64_t const cnt = 10;
	int64_t const iv[10][2] = {
		{ 100, 200 },
		{ 50, 100 },
		{ 120, 180 },
		{ 60, 240 },
		{ 70, 200 },
		{ 100, 140 },
		{ 90, 150 },
		{ 110, 160 },
		{ 120, 140 },
		{ 40, 190 }
	};

	/* insert */
	for(int64_t i = 0; i < cnt; i++) {
		struct ut_ivnode_s *n = (struct ut_ivnode_s *)
			ivtree_create_node(tree);
		n->h.lkey = iv[i][0];
		n->h.rkey = iv[i][1];
		n->val = 0x12345678;
		ivtree_insert(tree, (ivtree_node_t *)n);

		debug("print tree");
		ivtree_walk(tree, ivtree_print_node, NULL);
	}

	/* create iterator */
	ivtree_iter_t *iter = ivtree_contained(tree, 100, 170);
	assert(iter != NULL);

	ivtree_node_t *n = NULL;
	n = ivtree_next(iter); assert(_check_node(n, 100, 140), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 110, 160), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 120, 140), _print_node(n));
	n = ivtree_next(iter); assert(n == NULL);

	/* cleanup */
	ivtree_iter_clean(iter);
	ivtree_clean(tree);
}

/* containing iterator */
unittest()
{
	ivtree_t *tree = ivtree_init(sizeof(struct ut_ivnode_s), NULL);

	int64_t const cnt = 10;
	int64_t const iv[10][2] = {
		{ 100, 200 },
		{ 50, 100 },
		{ 120, 180 },
		{ 60, 240 },
		{ 70, 200 },
		{ 100, 140 },
		{ 90, 150 },
		{ 110, 160 },
		{ 120, 140 },
		{ 40, 190 }
	};

	/* insert */
	for(int64_t i = 0; i < cnt; i++) {
		struct ut_ivnode_s *n = (struct ut_ivnode_s *)
			ivtree_create_node(tree);
		n->h.lkey = iv[i][0];
		n->h.rkey = iv[i][1];
		n->val = 0x12345678;
		ivtree_insert(tree, (ivtree_node_t *)n);

		debug("print tree");
		ivtree_walk(tree, ivtree_print_node, NULL);
	}

	/* create iterator */
	ivtree_iter_t *iter = ivtree_containing(tree, 100, 150);
	assert(iter != NULL);

	ivtree_node_t *n = NULL;
	n = ivtree_next(iter); assert(_check_node(n, 40, 190), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 60, 240), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 70, 200), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 90, 150), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 100, 200), _print_node(n));
	n = ivtree_next(iter); assert(n == NULL);

	/* cleanup */
	ivtree_iter_clean(iter);
	ivtree_clean(tree);
}

/* intersect iterator */
unittest()
{
	ivtree_t *tree = ivtree_init(sizeof(struct ut_ivnode_s), NULL);

	int64_t const cnt = 10;
	int64_t const iv[10][2] = {
		{ 100, 200 },
		{ 50, 100 },
		{ 120, 180 },
		{ 60, 240 },
		{ 70, 200 },
		{ 100, 140 },
		{ 90, 150 },
		{ 110, 160 },
		{ 120, 140 },
		{ 40, 190 }
	};

	/* insert */
	for(int64_t i = 0; i < cnt; i++) {
		struct ut_ivnode_s *n = (struct ut_ivnode_s *)
			ivtree_create_node(tree);
		n->h.lkey = iv[i][0];
		n->h.rkey = iv[i][1];
		n->val = 0x12345678;
		ivtree_insert(tree, (ivtree_node_t *)n);

		debug("print tree");
		ivtree_walk(tree, ivtree_print_node, NULL);
	}

	/* create iterator */
	ivtree_iter_t *iter = ivtree_intersect(tree, 100, 120);
	assert(iter != NULL);

	ivtree_node_t *n = NULL;
	n = ivtree_next(iter); assert(_check_node(n, 40, 190), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 60, 240), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 70, 200), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 90, 150), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 100, 200), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 100, 140), _print_node(n));
	n = ivtree_next(iter); assert(_check_node(n, 110, 160), _print_node(n));
	n = ivtree_next(iter); assert(n == NULL);

	/* cleanup */
	ivtree_iter_clean(iter);
	ivtree_clean(tree);
}

/**
 * end of tree.c
 */
