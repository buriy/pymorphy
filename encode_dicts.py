#coding: utf-8
"""
Скрипт кодирования словарей. Пока только руками его править, если какие-то
настройки нужны.
"""

import codecs
import os

from pymorphy.backends import MrdDataSource, ShelveDataSource, PickleDataSource
from pymorphy.backends.trie_source import TrieDataSource

def _unlink(fn):
    try:
        os.unlink(fn)
    except OSError:
        pass

def convert_file(in_file, out_file, in_charset, out_charset):
    text = codecs.open(in_file, 'r', in_charset).read()
    codecs.open(out_file, 'w', out_charset ).write(text)


def load_mrd(dir):
    print('parsing source dictionary file...')
    mrd_source = MrdDataSource(os.path.join(dir,'morphs.utf8.mrd'),
                               os.path.join(dir,'gramtab.utf8.mrd'),
                               strip_EE=True)
    mrd_source.load()

    #print('calculating rule frequencies...')
    #mrd_source.calculate_rule_freq()
    return mrd_source


def make_pickled(dest_dir, mrd):
    print('creating pickled dictionary...')
    name = os.path.join(dest_dir, 'morphs.pickle')
    _unlink(name)
    source = PickleDataSource(dest_dir)
    source.convert_and_save(mrd)
    source.load()
    source._check_self()
    mrd._check_other(source)


def make_shelve(dest_dir, mrd, backend):
    print('creating %s dictionary...' % backend)

    for fn in ['lemmas', 'rules', 'endings', 'misc', 'freq']:
        _unlink(os.path.join(dest_dir, fn+'.'+backend))
        _unlink(os.path.join(dest_dir, fn+'.'+backend+'.tmp'))

    try:
        source = ShelveDataSource(dest_dir, backend)
        source.convert_and_save(mrd)
        source.load()
        mrd._check_other(source)
    except ImportError:
        print "Backend %s is not available." % backend


def make_trie(dest_dir, mrd, backend='trie', encoding='utf-8'):
    print('creating %s dictionary...' % backend)

    for fn in ['lemmas', 'rules', 'endings', 'misc', 'freq']:
        _unlink(os.path.join(dest_dir, fn+'.'+backend))
        _unlink(os.path.join(dest_dir, fn+'.'+backend+'.tmp'))

    try:
        source = TrieDataSource(dest_dir, backend, encoding=encoding)
        source.convert_and_save(mrd)
        source.load()
        #mrd._check_other(source)
    except ImportError:
        print "Backend %s is not available." % backend


def convert_dicts(src_dir, dest_dir, lang):

    if lang == 'en':
        msg = "encoding English.."
        src_subdir = os.path.join('SrcMorph','EngSrc')
        coding = 'latin1'
        gramtab = 'egramtab.tab'
    elif lang == 'ru':
        msg = "encoding Russian.."
        src_subdir = os.path.join('SrcMorph','RusSrc')
        coding = 'cp1251'
        gramtab = 'rgramtab.tab'
    else:
        print "invalid language"
        return

    print msg
    convert_file(os.path.join(src_dir, src_subdir, 'morphs.mrd'),
                 os.path.join(dest_dir, lang, 'morphs.utf8.mrd'),
                 coding, 'utf8')

    convert_file(os.path.join(src_dir, 'Morph', gramtab),
                 os.path.join(dest_dir, lang, 'gramtab.utf8.mrd'),
                 coding, 'utf8')


def cleanup_after_convert(dir):
    print('cleaning up...')
    os.unlink(os.path.join(dir, 'morphs.utf8.mrd'))
    os.unlink(os.path.join(dir, 'gramtab.utf8.mrd'))
    print("========")


if __name__ == '__main__':
    MrdDataSource.setup_psyco()

    src_dir = 'dicts/src/Dicts'
    dest_dir = 'dicts/converted'

# ======= en ===========
    convert_dicts(src_dir, dest_dir, 'en')
    en_dest = os.path.join(dest_dir, 'en')

    mrd = load_mrd(en_dest)
#    make_pickled(en_dest, mrd)
#    make_shelve(en_dest, mrd, 'shelve')      # bsddb
#    make_shelve(en_dest, mrd, 'cdb')         # cdb
#    make_shelve(en_dest, mrd, 'tch')         # tokyo cabinet hash
#    make_shelve(en_dest, mrd, 'tcb')         # tokyo cabinet btree+
#    make_shelve(en_dest, mrd, 'sqlite')      # sqlite
    make_trie(en_dest, mrd, encoding='latin1')# trie
    cleanup_after_convert(en_dest)

# ======= ru ===========
    convert_dicts(src_dir, dest_dir, 'ru')
    ru_dest = os.path.join(dest_dir, 'ru')

    mrd = load_mrd(ru_dest)
#    make_pickled(ru_dest, mrd)
#    make_shelve(ru_dest, mrd, 'shelve')     # bsddb
#    make_shelve(ru_dest, mrd, 'cdb')        # cdb
#    make_shelve(ru_dest, mrd, 'tch')        # tokyo cabinet hash
#    make_shelve(ru_dest, mrd, 'tcb')        # tokyo cabinet btree+
#    make_shelve(ru_dest, mrd, 'sqlite')     # sqlite
    make_trie(ru_dest, mrd, encoding='cp1251')# trie
    cleanup_after_convert(ru_dest)

    print "done."
