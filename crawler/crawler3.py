#-*- coding: utf-8 -*-
#import urllib.request as req
import requests as req
from urllib.parse import quote
import string
from bs4 import BeautifulSoup
import numpy as np
#https://blog.csdn.net/qq_25406563/article/details/81253347
#https://freelancerlife.info/zh/blog/python-web-scraping-user-agent-for-shopee/
#Request Headers->method: GET

#print(help(req))
def shopee_scraper(keyword,n_page=0,brands=False):
    #brands=4870#noCorrection=true
    #url = "https://shopee.tw/search?keyword="+keyword; #print(url)
    url = f'https://shopee.tw/search?keyword={keyword}&page={n_page}'
    if brands: url += '&brands=4870&noCorrection=true'
    url = quote(url, safe=string.printable).strip(); print(url)
    #Check for https://shopee.tw/robots.txt
    data = req.get(url, headers={
    'User-Agent':'Googlebot',
    'From':'jack661260@gmail.com',
    })
    print("HTTP status code: {}\n{}\n{}".format(data.status_code,data.history,data.url))
    root = BeautifulSoup(data.text, "html.parser"); #print(root)
    commodity = root.find_all("div",class_="_1w9jLI _37ge-4 _2ZYSiu")
    ID = np.zeros(shape=len(commodity),dtype=np.int32) 
    print("The number of commodities is {}".format(len(ID)))
    for i,com in enumerate(commodity):
        for ct in com.contents:
            if ct  == " - ": ID[i] = 1 
    all_items = root.find_all("div", class_="col-xs-2-4 shopee-search-item-result__item")
    ##links = [i.find('a').get('href') for i in all_items]
    links = [i.a.get('href') for i in all_items]
    contents = root.find_all("div",class_="_1NoI8_ _16BAGk")
    prices = root.find_all("span",class_="_341bF0")
    prices_range = []
    ii = 0
    for i,com in enumerate(commodity):
        if ID[i] == 0:
           prices_range.append(str(prices[ii].contents[0]))
           ii += 1      
        else:       
           prices_range.append(str(prices[ii].contents[0])+\
                         "-"+str(prices[ii+1].contents[0]))
           ii += 2      
    nn = 1
    for c,p,l in zip(contents, prices_range, links):
        print(str(nn)+"."+c.contents[0]); #print(c.string)
        print("$"+p)
        print('https://shopee.tw'+l)
        print('*---------------------------------*')
        nn += 1
#-------------------------------------------------------------------------  
keyword = "高蛋白"
shopee_scraper(keyword,1,True)
