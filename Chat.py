import streamlit as st
import os 
from Agent import Agent
from Tools import find_diseases_for_gene, downstream_analysis, collective_impact_analysys
from langchain_core.messages import AIMessage, HumanMessage
from config import OPENAI_API_KEY


assistant_img = 'assistant'
user_img = '👤'

with open('prompts/system_prompt.txt', 'r') as file:
    sys_prompt = file.read().replace('\n', '')

initial_message = ''' Hi! I'm an AI assistant built to understand and formulate hypotheses about gene involvement in
                        specific diseases. How can I help you today?'''


if "agent" not in st.session_state:
    st.session_state["agent"] = Agent(llm_provider='openai', 
                                llm_name='gpt-4o', 
                                api_key=OPENAI_API_KEY, 
                                temperature=0, 
                                tools=[find_diseases_for_gene, downstream_analysis, collective_impact_analysys], 
                                sys_prompt= sys_prompt,
                                initial_message=initial_message) 

agent = st.session_state["agent"]


for message in agent.chat_history:
    img = assistant_img if type(message) == AIMessage else user_img
    role = "assistant" if type(message) == AIMessage else "user"    
    with st.chat_message(role, avatar=img):
        st.markdown(message.content)


if prompt := st.chat_input(""):

    with st.chat_message("user", avatar=user_img):
        st.markdown(prompt)

    with st.chat_message("assistant", avatar=assistant_img):

        response = agent.chat(prompt)
        st.write(response)   
