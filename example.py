#coding: utf8

import cProfile, time
import re

import pymorphy

texts = []
texts.append(u'''
Сяпала Калуша по напушке и увазила бутявку. И волит:
— Калушата, калушаточки! Бутявка!
Калушата присяпали и бутявку стрямкали. И подудонились.
А Калуша волит:
— Оее, оее! Бутявка-то некузявая!
Калушата бутявку вычучили.
Бутявка вздребезнулась, сопритюкнулась и усяпала с напушки.
А Калуша волит:
— Бутявок не трямкают. Бутявки дюбые и зюмо-зюмо некузявые. От бутявок дудонятся.
А бутявка волит за напушкой:
— Калушата подудонились! Калушата подудонились! Зюмо некузявые! Пуськи бятые!
''')

texts.append(u'''
    Верьте аль не верьте,  а  жил  на  белом  свете  Федот-стрелец,  удалой
молодец.  Был  Федот ни красавец, ни урод, ни румян, ни бледен, ни богат, ни
беден, ни в парше, ни в парче, а так, вообче. Служба у Федота -- рыбалка  да
охота.  Царю  --  дичь  да  рыба, Федоту -- спасибо. Гостей во дворце -- как
семян в огурце. Один из Швеции, другой из Греции, третий с Гавай --  и  всем
жрать подавай! Одному -- омаров, другому -- кальмаров, третьему -- сардин, а
добытчик  один!  Как-то  раз  дают  ему  приказ: чуть свет поутру явиться ко
двору. Царь на вид сморчок, башка с кулачок, а злобности в ем --  агромадный
объем.  Смотрит  на  Федьку,  как  язвенник  на  редьку. На Федьке от страха
намокла рубаха, в висках застучало, в пузе заурчало, тут, как  говорится,  и
сказке начало...
    ''')

texts.append(u'''
Новые астрономические данные показали, что масса и скорость вращения
Млечного Пути оказались в несколько раз больше, чем считалось до сих пор.
Это известие, в свою очередь, означает, что наша Галактика столкнется со своей
соседкой - Туманностью Андромеды - раньше рассчитанного срока в пять
миллиардов лет. О своем открытии ученые доложили на встрече Американского
астрономического общества в Калифорнии.

В своих исследования специалисты из Гарвард-Смитсоновского астрофизического
центра и их коллеги использовали распределенную сеть радиотелескопов
(Very Long Baseline Array - VLBA). Они сосредоточились на наблюдении областей
интенсивного звездообразования. В некоторых районах внутри этих областей
находятся источники мазерного излучения. Электромагнитное излучение этого
типа возникает, когда облака межзвездного газа получают дополнительную энергию
от космических лучей. Получив такой "приток сил", молекулы газа переходят в
возбужденное состояние, а затем вновь возвращаются к "обычной жизни",
испуская при этом излишки энергии в виде фотонов.

Астрономы провели серию периодических наблюдений за выбранными двадцатью
источниками мазерного излучения. Такой подход позволил ученым определить
небольшие изменения положения мазеров относительно удаленных неподвижных звезд.
Все изучаемые мазеры находились на рукавах нашей спиральной Галактики, и
исследователи смогли определить ее структуру, построив трехмерную карту движения
источников мазерного излучения.

Полученные ими данные позволили оценить такой параметр Галактики, как скорость
вращения ее рукавов (в той области, где находится Солнечная система). Оказалось,
что Млечный Путь "крутится" со скоростью около 254 километров в секунду. Эта
цифра на 15% больше всех предыдущих значений. Скорость вращения галактик связана
с их массой, поэтому авторы работы смогли оценить и этот параметр Млечного Пути.
Согласно их расчетам, Галактика приблизительно в два раза тяжелее, чем считалось
ранее. А значит, силы притяжения между потяжелевшей до трех триллионов солнечных
масс Галактикой и Туманностью Андромеды также больше, чем принято считать. На
данный момент астрономы затрудняются назвать новую "дату" столкновения.

Еще одним открытием, которое было сделано с помощью VLBA, стало изменение число
рукавов Млечного Пути. Авторы исследования утверждают, что их не два, а четыре.

Данные, указывающие на то, что имеющаяся у астрономов информация о массе,
скорости и числе рукавов не полна, появляются не первый раз. Доказательства,
полученные в данной работе, являются более убедительными, чем предыдущие, так
как авторам удалось построить трехмерную карту структуры Галактики. Тем не менее,
как отмечает Роберт Бенджамин (Robert Benjamin) из Университета Висконсина
(University of Wisconsin-Whitewater), исследование ученых не ставит точку в вопросе
изучения нашего "космического дома", а всего лишь указывает направление дальнейшей работы.
    ''')

texts.append(u'''
Портрет кисти Фрэнсиса Бэкона под названием "Мужчина в синем VI" выставлен на
лондонские торги Christie's с оценкой 6 миллионов фунтов стерлингов (8,8 миллиона
долларов). Аукцион, сообщает The Times, назначен на 11 февраля 2009 года.
Работа датируется 1954 годом. С 1971-го она находилась в одних руках; нынешний
владелец купил ее за 31,5 тысячи фунтов.

Считается, что вопреки обыкновению Бэкон написал этот портрет с натуры (обычно
художник пользовался фотографией). "Мужчина в синем VI" входит в серию из семи
работ, которые Бэкон написал, поссорившись со своим любовником Питером Лэйси,
военным летчиком, и сбежав от него в гостиницу. Там художник и встретил некоего
постояльца, согласившегося позировать.

Живопись Фрэнсиса Бэкона (1909-1992) с 2006 года стремительно дорожала. В мае
2008-го его триптих ушел с молотка за 86,2 миллиона долларов, тремя месяцами
раньше "Обнаженная" была продана почти за 40 миллионов. Бэкон является одним
из самых востребованных и дорогих художников второй половины XX века.
    ''')
r = re.compile('[\W+-]',re.U)
words = []
for text in texts:
    words.extend(r.split(text.upper()))

pymorphy.setup_psyco()
morph = pymorphy.get_shelve_morph('ru')#, 'cdb', cached=True)
#morph = pymorphy.get_pickle_morph('ru')

print len(words)
def prof(words):
    for word in words:
        if word:
            norm_forms = morph.normalize(word)
            try:
                plural_forms = morph.pluralize_ru(word)
            except KeyError:
                print word
                raise
#            print forms

#for x in range(0,20):
#prof(words)

cProfile.run('prof(words)', sort='time')

time.sleep(10)
