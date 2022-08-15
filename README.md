# SIR 모델을 사용한 코로나19분석

4조 정지원, 임민석, 최대성, 최민준

SIR모델은 미분방정식을 이용해 감염병 확산을 파악하는데 사용하는 방법이다. 이때 SIR은 Susceptible Infected Recovered의 약자로 이를 이용한 다양한 발전된 모듈들이 존재한다.

이때 사용 되는 용어들은 아래와 같다

#R0은 기초 감염재생산수로 쉽게 말하면 최초감염자가 만들어내는 2차 감염자의 수로, R0이 1보다 크면 펜데믹, 1이라면 풍토병, 1보다 작으면 감염자가 줄어든다

#Y는 회복률로 감염기간의 역수로 정의 한다

#B는 감염의 효과율로 R0 * B로 정의 된다

이러한 모델에 잠복기를 추가한 SEIR모듈을 이용해 예측 했을때 2019년 12월 31일 기준으로 2025년 6월 22일에 종식할것이라고 예측 했다

다음으로 재감염준을 고려한 개량형 SEIRS모델을 사용 했을때 2019년 12월 31일 기준으로 2025년 6월 22일에 종식 할것으로 예측 했다

그러나 이는 어디 까지나 여러 요소들을 통한 예측이지만 실제 지표와는 큰차이를 보였다.
