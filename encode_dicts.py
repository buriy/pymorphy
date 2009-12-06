#coding: utf-8

import codecs
import os

from pymorphy.backends import MrdDataSource, ShelveDataSource, PickleDataSource

def convert_file(in_file, out_file, in_charset, out_charset):
    text = codecs.open(in_file, 'r', in_charset).read()
    codecs.open(out_file, 'w', out_charset ).write(text)

def create_shelf_and_pickle(lang):
    dir = 'dicts/converted/'+lang

    print('removing old data...')
    try:
        os.unlink(os.path.join(dir, 'morphs.pickle'))
    except OSError:
        pass

    try:
        os.unlink(os.path.join(dir, 'lemmas.shelve'))
        os.unlink(os.path.join(dir, 'rules.shelve'))
        os.unlink(os.path.join(dir, 'endings.shelve'))
        os.unlink(os.path.join(dir, 'misc.shelve'))
        os.unlink(os.path.join(dir, 'freq.shelve'))
    except OSError:
        pass

    print('parsing source dictionary file...')
    mrd_source = MrdDataSource(os.path.join(dir,'morphs.utf8.mrd'),
                       os.path.join(dir,'gramtab.utf8.mrd'),
                       strip_EE=True)
    mrd_source.load()

    print('calculating rule frequencies...')
    mrd_source.calculate_rule_freq()

    print('creating pickled dictionary...')
    pickle_source = PickleDataSource(os.path.join(dir,'morphs.pickle'))
    pickle_source.convert_and_save(mrd_source)

    print('creating shelve dictionary...')
    shelve_source = ShelveDataSource(dir)
    shelve_source.convert_and_save(mrd_source)

    print('cleaning up...')
    os.unlink(os.path.join(dir,'morphs.utf8.mrd'))
    os.unlink(os.path.join(dir,'gramtab.utf8.mrd'))

    print "========"

if __name__ == '__main__':

    src_dir = 'dicts/src/Dicts'
    dest_dir = 'dicts/converted'

    print "encoding english.."
    convert_file(os.path.join(src_dir,'SrcMorph/EngSrc/morphs.mrd'),
                 os.path.join(dest_dir,'en/morphs.utf8.mrd'), 'latin1', 'utf8')
    convert_file(os.path.join(src_dir,'Morph/egramtab.tab'),
                 os.path.join(dest_dir,'en/gramtab.utf8.mrd'), 'latin1', 'utf8')

    create_shelf_and_pickle('en')

    print "encoding russian.."
    convert_file(os.path.join(src_dir,'SrcMorph/RusSrc/morphs.mrd'),
                 os.path.join(dest_dir,'ru/morphs.utf8.mrd'), 'cp1251', 'utf8')
    convert_file(os.path.join(src_dir,'Morph/rgramtab.tab'),
                 os.path.join(dest_dir,'ru/gramtab.utf8.mrd'), 'cp1251', 'utf8')

    create_shelf_and_pickle('ru')

    print "done."
