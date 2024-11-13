import streamlit as st

#1:ファイルの読み込み
from Bio import SeqIO

#ファイルの読み込み 
record1 = next(SeqIO.parse("NC_045512_Wuhan.fa", "fasta"))
record2 = next(SeqIO.parse("sars.fa", "fasta"))

#配列の取り出し
seq1 = record1.seq
seq2 = record2.seq

#3:画像の表示
import numpy as np
import matplotlib.pyplot as plt

win = 10
len1 = len(seq1) - win + 1 #配列１の長さ
len2 = len(seq2) - win + 1 #配列２の長さ

width = 500 #実際に描く幅
height = 500 #実際に描く高さ

image= np. zeros((height, width)) #実際に描く幅・高さの行列

#ハッシュの生成
hash = {}

for x in range(len1):
    sub1 = seq1[x : x + win] #seq1からwinまでの部分配列を作成
    if sub1 not in hash: #hashにsub1が存在しないとき
        hash[sub1]=[] #空のリストを登録
    hash[sub1].append(x) #位置xを追加

for y in range(len2):
    sub2 = seq2[y : y + win] #seq2からwinまでの部分配列を表示
    py = int(y/len2*height) #yを画像位置pyに
    if sub2 in hash:
        for x in hash[sub2]: #sub2に対応するseq1内のすべての位置xに対して
            px = int(x/len1*width) #xを画像位置pxに
            image[py, px] = 1 #[py, px]にドット


plt. imshow( image, extent=(1,len1,len2,1),cmap="Grays")
#カラーマップを灰色濃淡
#0としかないので白黒になる
#plt. show()

# plt.show() 最 後の行をコメントアウト

st.pyplot(plt) # 書 き 加える：Streamlit上にMatplotlibを表示