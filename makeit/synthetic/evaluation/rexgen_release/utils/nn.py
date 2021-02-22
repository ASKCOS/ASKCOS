import numpy as np
from scipy.sparse import csc_matrix, vstack
import tensorflow.compat.v1 as tf
import math

def batch_normalization(x, scope, decay=0.999, eps=1e-6, training=True):
    ndim = len(x.get_shape().as_list())
    fdim = x.get_shape().as_list()[-1]
    with tf.variable_scope(scope):
        gamma = tf.get_variable("scale", [fdim], tf.float32, tf.constant_initializer(1.0))
        beta = tf.get_variable("offset", [fdim], tf.float32, tf.constant_initializer(0.0))
        mean = tf.get_variable("mean", [fdim], tf.float32, tf.constant_initializer(0.0), trainable=False)
        var = tf.get_variable("variance", [fdim], tf.float32, tf.constant_initializer(1.0), trainable=False)
        if training:
            x_mean, x_var = tf.nn.moments(x, list(range(ndim - 1)))
            avg_mean = tf.assign(mean, mean * decay + x_mean * (1.0 - decay))
            avg_var = tf.assign(var, var * decay + x_var * (1.0 - decay))
            with tf.control_dependencies([avg_mean, avg_var]):
                return tf.nn.batch_normalization(x, x_mean, x_var, beta, gamma, eps)
        else:
            return tf.nn.batch_normalization(x, mean, var, beta, gamma, eps)

def batch_normalization_with_mask(x, mask, scope, decay=0.999, eps=1e-6, training=True):
    ndim = len(x.get_shape().as_list())
    fdim = x.get_shape().as_list()[-1]
    with tf.variable_scope(scope):
        gamma = tf.get_variable("scale", [fdim], tf.float32, tf.constant_initializer(1.0))
        beta = tf.get_variable("offset", [fdim], tf.float32, tf.constant_initializer(0.0))
        mean = tf.get_variable("mean", [fdim], tf.float32, tf.constant_initializer(0.0), trainable=False)
        var = tf.get_variable("variance", [fdim], tf.float32, tf.constant_initializer(1.0), trainable=False)
        if training:
            x_mean, x_var = tf.nn.weighted_moments(x, list(range(ndim - 1)), mask)
            avg_mean = tf.assign(mean, mean * decay + x_mean * (1.0 - decay))
            avg_var = tf.assign(var, var * decay + x_var * (1.0 - decay))
            with tf.control_dependencies([avg_mean, avg_var]):
                return tf.nn.batch_normalization(x, x_mean, x_var, beta, gamma, eps)
        else:
            return tf.nn.batch_normalization(x, mean, var, beta, gamma, eps)

def smart_gatherND(params, indices):
    shape = indices.get_shape()
    ndim = len(shape)
    cur_shape = tf.shape(indices)
    new_shape = tf.gather(cur_shape, list(range(ndim)))
    new_shape = tf.concat(0, [new_shape, hidden_size])
    indices = tf.reshape([-1])
    res = tf.gather(params, indices)
    res = tf.reshape(new_shape)
    res.set_shape(shape[:-1] + [hidden_size])
    return res

def lookup_table(input_, vocab_size, output_size, scope):
    with tf.variable_scope(scope):
        embedding = tf.get_variable("embedding", [vocab_size, output_size], tf.float32, tf.random_normal_initializer(stddev=1.0 / math.sqrt(output_size)))
    return tf.nn.embedding_lookup(embedding, input_)

def sparse_linear(input_, input_size, output_size, scope, init_bias=0.0):
    with tf.variable_scope(scope):
        W = tf.get_variable("Matrix", [input_size, output_size], tf.float32, tf.random_normal_initializer(stddev=1.0 / math.sqrt(output_size)))
        b = tf.get_variable("bias", [output_size], initializer=tf.constant_initializer(init_bias))
    return tf.sparse_tensor_dense_matmul(input_, W) + b

def linear(input_, output_size, scope, reuse=False, init_bias=0.0):
    shape = input_.get_shape().as_list()
    stddev = min(1.0 / math.sqrt(shape[-1]), 0.1)
    with tf.variable_scope(scope, reuse=reuse):
        W = tf.get_variable("Matrix", [shape[-1], output_size], tf.float32, tf.random_normal_initializer(stddev=stddev))
    if init_bias is None:
        return tf.matmul(input_, W)
    with tf.variable_scope(scope, reuse=reuse):
        b = tf.get_variable("bias", [output_size], initializer=tf.constant_initializer(init_bias))
    return tf.matmul(input_, W) + b

def linearND(input_, output_size, scope, reuse=False, init_bias=0.0):
    shape = input_.get_shape().as_list()
    ndim = len(shape)
    stddev = min(1.0 / math.sqrt(shape[-1]), 0.1)
    with tf.variable_scope(scope, reuse=reuse):
        W = tf.get_variable("Matrix", [shape[-1], output_size], tf.float32, tf.random_normal_initializer(stddev=stddev))
    X_shape = tf.gather(tf.shape(input_), list(range(ndim-1)))
    target_shape = tf.concat([X_shape, [output_size]], 0)
    exp_input = tf.reshape(input_, [-1, shape[-1]])
    if init_bias is None:
        res = tf.matmul(exp_input, W)
    else:
        with tf.variable_scope(scope, reuse=reuse):
            b = tf.get_variable("bias", [output_size], initializer=tf.constant_initializer(init_bias))
        res = tf.matmul(exp_input, W) + b
    res = tf.reshape(res, target_shape)
    res.set_shape(shape[:-1] + [output_size])
    return res

def CSR2TF(x):
    indices = [[0,0]] # Avoid empty matrix
    values = [0.0]
    for i in range(x.shape[0]):
        left = x.indptr[i]
        right = x.indptr[i + 1]
        for j in range(left, right):
            indices.append([i, x.indices[j]])
            values.append(x.data[j])
    return (indices, values, x.shape)

