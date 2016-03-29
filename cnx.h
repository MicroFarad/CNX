/*
Header file for signal processing algorithm connector data structures
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#ifndef CNX_H
#define CNX_H

#ifndef CNX_Scalar
#define CNX_Scalar float
#endif

#define CNX_THREADING

#ifdef CNX_THREADING
#include <threads.h>
#define CNX_THREADED   1
#define CNX_POOLBOUND  3
#define CNX_QUEUEBOUND 5
#endif

#ifndef NULL
#define NULL ((void*)0)
#endif

typedef int (*CNX_Algorithm)(void*, int);

typedef struct CNX_Buffer
{
	void *data;
	struct CNX_Buffer *next;
} CNX_Buffer;

typedef struct
{
	CNX_Buffer *pool;
	CNX_Buffer *queueold;
	CNX_Buffer *queuenew;
	int bufferwidth;
	int bufferlength;
#ifdef CNX_THREADING
	int flags;
	mtx_t poolmutex;
	mtx_t queuemutex;
	cnd_t poolcond;
	cnd_t queuecond;
#endif
} CNX_Connector;

#ifdef CNX_THREADING
typedef struct
{
	CNX_Algorithm algorithm;
	thrd_t thread;
	void *data;
	int go;
} CNX_Instance;

CNX_Instance *CNX_StartInstance(CNX_Algorithm algorithm, void *data);

int CNX_StopInstance(CNX_Instance *instance);
#endif

void CNX_ConstructConnector(CNX_Connector *connector, int bufferwidth, int bufferlength, int poolsize, int flags);

CNX_Buffer *CNX_TakePool(CNX_Connector *connector);

void CNX_QueueBuffer(CNX_Connector *connector, CNX_Buffer *buffer);

CNX_Buffer *CNX_TakeQueue(CNX_Connector *connector);

void CNX_PoolBuffer(CNX_Connector *connector, CNX_Buffer *buffer);

void CNX_CleanPool(CNX_Connector *connector);

void CNX_DestroyConnector(CNX_Connector *connector);

#endif
