#include <stdlib.h>
#include <jni.h>

char *GetStringNativeChars(JNIEnv *env, jstring jstr);

JNIEXPORT jint JNICALL Java_ffe_core_FFEExec_NativeExec
   (JNIEnv *env, jobject obj, jstring jargv,
                              jstring jdir,
                              jstring jpath,
                              jstring jclasspath,
                              jstring jld_library_path) {
   const jbyte *argv;
   const jbyte *dir;
   const jbyte *path;
   const jbyte *classpath;
   char *temp;
   char *current;
   int len;

   // Convert Java Strings to jbyte arrays
   path = GetStringNativeChars(env, jpath);
   classpath = GetStringNativeChars(env, jclasspath);
   dir = GetStringNativeChars(env, jdir);
   argv = GetStringNativeChars(env, jargv);

   // Set the PATH

   current = getenv("PATH");
   if (current != NULL) {
      temp = strstr(current,path);
      if (temp == NULL) {
         len  = strlen(path) + strlen(current) + 20;
         temp = (char *) malloc(len);
         sprintf(temp, "PATH=%s;%s", path, current);
         putenv(temp);
         free(temp);
      }
      free(current);
   } else {
      temp = (char *) malloc(strlen(path) + 20);
      sprintf(temp, "PATH=%s", path);
      putenv(temp);
      free(temp);
   }

   // Set the CLASSPATH
   current = getenv("CLASSPATH");
   if (current != NULL) {
      temp = strstr(current,classpath);
      if (temp == NULL) {
         len  = strlen(classpath) + strlen(current) + 20;
         temp = (char *) malloc(len);
         sprintf(temp, "CLASSPATH=%s;%s", classpath, current);
         putenv(temp);
         free(temp);
      }
      free(current);
   } else {
      temp = (char *) malloc(strlen(classpath) + 20);
      sprintf(temp, "CLASSPATH=%s", classpath);
      putenv(temp);
      free(temp);
   }

   //printf("\n%s\n\n%s\n", getenv("PATH"), getenv("CLASSPATH"));
   //printf("%s\n", argv);

   chdir(dir);
   return system(argv);
}

char *GetStringNativeChars(JNIEnv *env, jstring jstr) {
   jbyteArray bytes = 0;
   jthrowable exc;
   jclass class;
   static jmethodID MID_String_getBytes = NULL;
   char *result = 0;

   if ((*env)->EnsureLocalCapacity(env, 2) < 0) {
      printf("Out of Memory\n");
      return 0; /* Out of Memory */
   }

   /* Do we have the cached Method ID Yet? */
   if (MID_String_getBytes == NULL) {
      class = (*env)->FindClass(env, "Ljava/lang/String;");
      if (class == 0) {
         printf("Could not find String Class\n");
         return 0; /* Could not Find the String Class */
      }
      MID_String_getBytes = (*env)->GetMethodID(env, class, "getBytes", "()[B");
      if (MID_String_getBytes == 0) {
         printf("Could not find getBytes\n");
         return 0;
       }
   }

   bytes = (*env)->CallObjectMethod(env, jstr, MID_String_getBytes);
   exc = (*env)->ExceptionOccurred(env);
   if (!exc) {
      jint len = (*env)->GetArrayLength(env, bytes);
      result = (char *) malloc(len + 1);
      if (result == 0) {
         (*env)->DeleteLocalRef(env, bytes);
         return 0;
      }
      (*env)->GetByteArrayRegion(env,bytes,0,len,(jbyte *)result);
      result[len] = 0;
   } else {
      printf("Exception calling String->GetBytes\n");
      (*env)->DeleteLocalRef(env, exc);
   }
   (*env)->DeleteLocalRef(env, bytes);
   return result;
}
