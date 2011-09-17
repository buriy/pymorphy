#coding: utf-8
from chartrie import CharTrie, FrozenCharTrie
from pymorphy.backends.base import DictDataSource
import os
from pymorphy.backends.shelve_source.shelf_with_hooks import ShelfWithHooks

DEBUG = True

class TrieDataSource(DictDataSource):
    """ Источник данных для морфологического анализатора pymorphy,
        берущий информацию из key-value базы данных, используя модифицированный
        интерфейс shelve из стандартной библиотеки. Позволяет не держать все
        данные в памяти и в то же время обеспечивает достаточно быструю скорость
        работы.

        Грамматическая информация и префиксы загружаются в память сразу.
    """

    def __init__(self, path='', db_type=None, encoding='utf-8'):
        self.path = path
        self.db_type = db_type or 'shelve'
        self.cached = False
        self.encoding = encoding

        super(TrieDataSource, self).__init__()

    def _load_file(self, name):
        in_ = open(self._path(name), 'rb')
        contents = in_.read()
        in_.close()
        return contents

    def _save_file(self, name, contents):
        out_ = open(self._path(name), 'wb')
        out_.write(contents)
        out_.close()

    def load(self):
        print "Using trie backend (beta)"
        print "Loading lemmas"
        self.lemmas_trie = FrozenCharTrie()
        self.lemmas_trie.loads(self._load_file('lemmas')) 
        print "Loading suffixes"
        self.suffixes_trie = FrozenCharTrie()
        self.suffixes_trie.loads(self._load_file('suffixes'))
        print "Loading endings"
        
        self.rules   = self._get_shelf('rules', 'r', 'int', 'pickle')
        self.endings = self._get_shelf('endings', 'r', 'unicode')

        misc = self._get_shelf('misc', 'r', 'unicode', 'pickle')
        self.gramtab = misc['gramtab']

        self.prefixes = set(misc['prefixes'])
        self.possible_rule_prefixes = set(misc['possible_rule_prefixes'])

        self.paradigm_list = misc['paradigm_list']
        self.suffix_dict = misc['suffix_dict']
        self.encoding = misc['encoding']
        self.initial_forms = misc['initial_forms']
        
    def convert_and_save(self, data_obj):
        print "Filling lemmas tree..."
        enc = self.encoding
        lemmas = CharTrie()
        paradigm_list = []
        paradigm_dict = {}
        for lemma, paradigms in data_obj.lemmas.iteritems():
            if lemma == '#':
                lemma = ''
            paradigms.sort()
            tupar = tuple(paradigms)
            try:
                lemma_id = paradigm_dict[tupar]
            except KeyError:
                paradigm_list.append(tupar)  
                lemma_id = paradigm_dict[tupar] = len(paradigm_list) - 1 
                
            lemmas[lemma.encode(enc)] = lemma_id
        if DEBUG:
            print "Lemma trie size:", len(lemmas)
            print "Paradigm count:", len(paradigm_list)

        print "Filling suffixes tree..."
        from collections import defaultdict
        suffix_dict = defaultdict(list)
        suffix_count = 0
        suffixes = CharTrie()
        
        initial_forms = {} 
        for paradigm_id, rules in data_obj.rules.iteritems():
            initial_forms[paradigm_id] = rules[0][0]
            for rule_suffix, rule_ancode, rule_prefix in rules:
                enc_suffix = rule_suffix.encode(enc)[::-1]
                suffix_id = suffixes[enc_suffix]
                if suffix_id is None:
                    suffix_id = suffixes[enc_suffix] = suffix_count
                    suffix_count += 1  
                suffix_dict[suffix_id, paradigm_id].append((rule_ancode, rule_prefix))

        if DEBUG:
            print "Done."
            para_count = sum(map(len, data_obj.lemmas.itervalues()))
            rule_count = sum(map(len, data_obj.rules.itervalues()))
            rule_count2 = sum(map(len, suffix_dict.itervalues()))
            assert rule_count == rule_count2
            print "Total number of paradigms: %s" % para_count
            print "Total number of rules: %s" % rule_count
            print "Total number of rules2: %s" % rule_count2
    
            print "Suffix trie size:", len(suffixes)
            print "Comb size:", len(suffix_dict)

            if False: # ultra-hardcore-debug mode
                from pymorphy.console import reprint
                
                for k,v in data_obj.__dict__.iteritems():
                    if isinstance(v, dict) and len(v)>=3:
                        print "Length of %s = %s" % (k, len(v)),
                        k2 = v.keys()[0]
                        print "Sample:", 
                        reprint((k2, '->', v[k2]))
                    elif isinstance(v, (list, tuple, set)) and len(v)>=3:
                        print "Length of %s = %s" % (k, len(v)),
                        print "Sample 0:", 
                        reprint(list(v)[0])
                    else:
                        print "Value of  %s =" % k,
                        reprint(v)

        self._save_file('lemmas', lemmas.dumps())
        self._save_file('suffixes', suffixes.dumps())

        endings_shelve = self._get_shelf('endings', 'c', 'unicode')

        for end, value in data_obj.endings.iteritems():
            endings_shelve[end] = value
        endings_shelve.close()

        rules_shelve = self._get_shelf('rules', 'c', 'int', 'pickle')
        for rule in data_obj.rules:
            rules_shelve[rule] = data_obj.rules[rule]
        rules_shelve.close()

        misc_shelve = self._get_shelf('misc', 'c', 'unicode', 'pickle')
        misc_shelve['encoding'] = enc
        misc_shelve['gramtab'] = data_obj.gramtab
        misc_shelve['prefixes'] = list(data_obj.prefixes)
        misc_shelve['possible_rule_prefixes'] = list(data_obj.possible_rule_prefixes)
        misc_shelve['paradigm_list'] = paradigm_list
        misc_shelve['suffix_dict'] = suffix_dict
        misc_shelve['initial_forms'] = initial_forms
        misc_shelve.close()

#        if data_obj.rule_freq:
#            freq_shelve = self._get_shelf('freq', 'c', 'int')
#            for (rule, freq,) in data_obj.rule_freq.items():
#                freq_shelve[int(rule)] = freq
#            freq_shelve.close()
#

    def _path(self, name):
        return os.path.join(self.path, name+'.'+self.db_type)

    def _get_shelf_class(self):
        if self.db_type == 'cdb':
            from cdb_shelve import CdbShelf
            return CdbShelf
        elif self.db_type == 'tch':
            from pytc_shelve import PytcHashShelf
            return PytcHashShelf
        elif self.db_type == 'tcb':
            from pytc_shelve import PytcBtreeShelf
            return PytcBtreeShelf
        elif self.db_type == 'sqlite':
            from sqlite_shelve import SqliteShelf
            return SqliteShelf
        return ShelfWithHooks

    def _get_shelf(self, filename, *args, **kwargs):
        path = self._path(filename)
        return self._get_shelf_class()(path, cached=self.cached,  *args, **kwargs)

    def _check_self(self):
        raise NotImplementedError()

    def _check_other(self, data_source):
        raise NotImplementedError()
