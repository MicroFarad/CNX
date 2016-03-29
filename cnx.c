/*
Source file for signal processing algorithm connector data structures
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#include "cnx.h"
#include <stdlib.h>

#ifdef CNX_THREADING
int RunContinuously(void *instance)
{
	int count = 0;
	int *go = &((CNX_Instance*)instance)->go;
	CNX_Algorithm algorithm = ((CNX_Instance*)instance)->algorithm;
	void *data = ((CNX_Instance*)instance)->data;
	while(*go) count += algorithm(data, 1);
	return count;
}

CNX_Instance *CNX_StartInstance(CNX_Algorithm algorithm, void *data)
{
	CNX_Instance *instance = malloc(sizeof(CNX_Instance));
	instance->algorithm = algorithm;
	instance->data = data;
	instance->go = 1;
	thrd_create(&instance->thread, RunContinuously, instance);
	return instance;
}

int CNX_StopInstance(CNX_Instance *instance)
{
	int count;
	instance->go = 0;
	thrd_join(instance->thread, &count);
	free(instance);
	return count;
}
#endif

void CNX_ConstructConnector(CNX_Connector *connector, int bufferwidth, int bufferlength, int poolsize, int flags)
{
	connector->pool = NULL;
	connector->queueold = NULL;
	connector->queuenew = NULL;
	connector->bufferwidth = bufferwidth;
	connector->bufferlength = bufferlength;
	for(int n = 0; n < poolsize; n++)
	{
		CNX_Buffer *buffer = malloc(sizeof(CNX_Buffer));
		buffer->data = malloc(bufferwidth*bufferlength);
		buffer->next = connector->pool;
		connector->pool = buffer;
	}
#ifdef CNX_THREADING
	connector->flags = flags;
	if(connector->flags & CNX_THREADED)
	{
		mtx_init(&connector->poolmutex, mtx_plain);
		mtx_init(&connector->queuemutex, mtx_plain);
	}
	if(connector->flags & CNX_POOLBOUND & ~CNX_THREADED) cnd_init(&connector->poolcond);
	if(connector->flags & CNX_QUEUEBOUND & ~CNX_THREADED) cnd_init(&connector->queuecond);
#endif
}

CNX_Buffer *CNX_TakePool(CNX_Connector *connector)
{
	CNX_Buffer *buffer;
#ifdef CNX_THREADING
	if(connector->flags & CNX_THREADED) mtx_lock(&connector->poolmutex);
	if(connector->flags & CNX_POOLBOUND & ~CNX_THREADED) while(!connector->pool) cnd_wait(&connector->poolcond, &connector->poolmutex);
#endif
	if(connector->pool)
	{
		buffer = connector->pool;
		connector->pool = buffer->next;
	}
	else
	{
		buffer = malloc(sizeof(CNX_Buffer));
		buffer->data = malloc(connector->bufferwidth*connector->bufferlength);
	}
#ifdef CNX_THREADING
	if(connector->flags & CNX_THREADED) mtx_unlock(&connector->poolmutex);
#endif
	return buffer;
}

void CNX_PoolBuffer(CNX_Connector *connector, CNX_Buffer *buffer)
{
#ifdef CNX_THREADING
	if(connector->flags & CNX_THREADED) mtx_lock(&connector->poolmutex);
#endif
	buffer->next = connector->pool;
	connector->pool = buffer;
#ifdef CNX_THREADING
	if(connector->flags & CNX_POOLBOUND & ~CNX_THREADED) cnd_signal(&connector->poolcond);
	if(connector->flags & CNX_THREADED) mtx_unlock(&connector->poolmutex);
#endif
}

CNX_Buffer *CNX_TakeQueue(CNX_Connector *connector)
{
	CNX_Buffer *buffer = NULL;
#ifdef CNX_THREADING
	if(connector->flags & CNX_THREADED) mtx_lock(&connector->queuemutex);
	if(connector->flags & CNX_QUEUEBOUND & ~CNX_THREADED) while(!connector->queueold) cnd_wait(&connector->queuecond, &connector->queuemutex);
#endif
	buffer = connector->queueold;
	if(buffer) connector->queueold = buffer->next;
	if(!buffer->next) connector->queuenew = NULL;
#ifdef CNX_THREADING
	if(connector->flags & CNX_THREADED) mtx_unlock(&connector->queuemutex);
#endif
	return buffer;
}

void CNX_QueueBuffer(CNX_Connector *connector, CNX_Buffer *buffer)
{
#ifdef CNX_THREADING
	if(connector->flags & CNX_THREADED) mtx_lock(&connector->queuemutex);
#endif
	if(connector->queuenew) connector->queuenew->next = buffer;
	else connector->queueold = buffer;
	connector->queuenew = buffer;
	buffer->next = NULL;
#ifdef CNX_THREADING
	if(connector->flags & CNX_QUEUEBOUND & ~CNX_THREADED) cnd_signal(&connector->queuecond);
	if(connector->flags & CNX_THREADED) mtx_unlock(&connector->queuemutex);
#endif
}

void CNX_CleanPool(CNX_Connector *connector)
{
#ifdef CNX_THREADING
	if(connector->flags & CNX_THREADED) mtx_lock(&connector->poolmutex);
#endif
	CNX_Buffer *buffer = connector->pool;
	while(buffer)
	{
		CNX_Buffer *next = buffer->next;
		free(buffer->data);
		free(buffer);
		buffer = next;
	}
	connector->pool = NULL;
#ifdef CNX_THREADING
	if(connector->flags & CNX_THREADED) mtx_unlock(&connector->poolmutex);
#endif
}

void CNX_DestroyConnector(CNX_Connector *connector)
{
	while(connector->queuenew) CNX_PoolBuffer(connector, CNX_TakeQueue(connector));
	CNX_CleanPool(connector);
#ifdef CNX_THREADING
	if(connector->flags & CNX_THREADED)
	{
		mtx_destroy(&connector->poolmutex);
		mtx_destroy(&connector->queuemutex);
	}
	if(connector->flags & CNX_POOLBOUND & ~CNX_THREADED) cnd_destroy(&connector->poolcond);
	if(connector->flags & CNX_QUEUEBOUND & ~CNX_THREADED) cnd_destroy(&connector->queuecond);
#endif
}
